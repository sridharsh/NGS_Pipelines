#!/usr/bin/env nextflow

ref_genome = file("reference_genome/GRCh38.fa")

params.reads = "/R{1,2}.fastq.gz"
read_files = file(params.reads)

/* ALIGNMENT using BWA for the fastq files (reads) and gives a bam file of the resulting alignment */

process aligning{
	
	input:
	ref_genome
	read_files

	output:
	file "Aligned.bam" into aligned_results

	shell:
	"""
	module load bwa/0.7.17 
	module load samtools/1.9
	bwa mem ${ref_genome} ${read_files[0]} ${read_files[1]} | samtools view -S -b -o Aligned.bam
	"""

}

process filter_sort{

	input:
	file "Aligned.bam" from aligned_results

	output:
	set "Aligned_Sorted.bam", "Aligned_Sorted_Filtered.bam", "Aligned_Sorted_Filtered.bai" into sam_results

	shell:
	'''
	module load samtools/1.9
	samtools sort Aligned.bam -o Aligned_Sorted.bam
	samtools view -b -F 1804 -q 2 Aligned_Sorted.bam -o Aligned_Sorted_Filtered.bam
	samtools index Aligned_Sorted_Filtered.bam
	'''
	
}

process Variant_Calling{

	input:
	set "Aligned.sorted.bam", "Aligned.sorted.filtered.bam" from sam_results

	output:
	file "Aligned.sorted.filtered.RG.bam"
	set "Aligned-SF-RG-Marked-Dups.bam", "Marked-Dup-Metrics.txt" into Marked_dups
	file "Raw_Variants.vcf" into Variants

	shell:
	"""
	module load picard/2.18.4
	java -jar $PICARD AddOrReplaceReadGroups I=Aligned.sorted.filtered.bam O=Aligned.sorted.filtered.RG.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
	java -jar $PICARD MarkDuplicates I=Aligned.sorted.filtered.RG.bam O=Aligned-SF-RG-Marked-Dups.bam M=Marked-Dup-Metrics.txt
	module load samtools/1.9
	samtools index Aligned-SF-RG-Marked-Dups.bam
	module unload java
	module load gatk/4.0.1.2
	gatk HaplotypeCaller -I Aligned-SF-RG-Marked-Dups.bam -ERC GVCF -O Raw_Variants.vcf -R ${ref_genome}
	"""
}
