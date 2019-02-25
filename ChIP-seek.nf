#!/usr/bin/env nextflow

ref_genome = file("GRCh38.fa")

params.reads = "reads*.fastq.gz"
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

process peak_calling{

        input:
        set "Aligned_Sorted.bam", "Aligned_Sorted_Filtered.bam", "Aligned_Sorted_Filtered.bai" from sam_results

        shell:
        '''
        module load macs/2.1.0
        macs2 callpeak -t Aligned_Sorted_Filtered.bam -f BAM -g hs -n test_nextflow -B -q 0.01 
        '''

}
