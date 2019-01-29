#!/usr/bin/env nextflow

ref_genome = file("hg38_Star_directory")
gtf_file = file("GRCh38_reference.gtf")
GrCh38 = file("GRCh38_reference_genome.fa")

params.reads = "./R{1,2}.fastq"
read_files = file(params.reads)

/* ALIGNMENT using STAR for the fastq files (reads) and gives a bam file of the resulting alignment */

process aligning{
	
	input:
	ref_genome
	read_files

	output:
	set "Aligned.out.sam", "Aligned.out.bam" into aligned_results

	shell:
	"""
	module unload gcc
	module load star/2.6.1d
	STAR --genomeDir ${ref_genome} --readFilesIn ${read_files[0]} ${read_files[1]} --runThreadN 10
	module load samtools/1.9
        samtools view Aligned.out.sam -b -S -o Aligned.out.bam
	"""
}

process filter_sort{

        input:
        set "Aligned.out.sam", "Aligned.out.bam" from aligned_results

        output:
        set "Aligned.sorted.bam", "Aligned.sorted.filtered.bam" into sam_results

        shell:
	"""
	module load samtools/1.9
	samtools sort Aligned.out.bam -o Aligned.sorted.bam
	samtools  view -b -O BAM -F 1804 -q 2 Aligned.sorted.bam -o Aligned.sorted.filtered.bam
	"""
}

process cufflinking{

	input:
	set "Aligned.sorted.bam", "Aligned.sorted.filtered.bam" from sam_results

	output:
	set "skipped.gtf", "isoforms.fpkm_tracking", "genes.fpkm_tracking", 'transcripts.gtf' into transcripts

	shell:
	"""
	module load cufflinks/2.2.1
	cufflinks Aligned.sorted.filtered.bam
	"""
}

process feature_counting{
	
	input:
	set "Aligned.sorted.bam", "Aligned.sorted.filtered.bam" from sam_results
	gtf_file	

	output:
	set "counts.txt", "counts.txt.summary" into feature_counts	

	shell:
	"""
	module load subread/1.6.3
	featureCounts Aligned.sorted.filtered.bam -o counts.txt -a ${gtf_file}
	"""
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
	gatk HaplotypeCaller -I Aligned-SF-RG-Marked-Dups.bam -ERC GVCF -O Raw_Variants.vcf -R ${GrCh38} 
	"""
}
