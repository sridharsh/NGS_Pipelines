#!/usr/bin/env nextflow

ref_genome = file("/sc/orga/projects/PBG/REFERENCES/GRCh38/BWAindex/GRCh38.primary_assembly.genome.fa")
params.reads = "/sc/orga/projects/PBG/hara/hara-tests/dnaseq/4680_6859_ACAGTG_L002_R{1,2}_002.HB676ADXX.fastq.gz"
read_files = file(params.reads)


/* ALIGNMENT using BWA for the fastq files (reads) and gives a bam file of the resulting alignment */

process aligning{

        executor = 'lsf'
        queue = 'premium'
        clusterOptions = '-P acc_apollo -sp 50'
        cpus = 1
        memory = "12GB" // 4*12
        time = "40h"

        input:
        ref_genome
        read_files

        output:
        file "Aligned.bam" into aligned_bam

        module "bwa/0.7.17:samtools/1.9"

        shell:
        """
        bwa mem ${ref_genome} ${read_files[0]} ${read_files[1]} | samtools view -S -b -o Aligned.bam
        """

}

process filter_sort{

        executor = 'lsf'
        queue = 'premium'
        clusterOptions = '-P acc_apollo -sp 50'
        cpus = 1
        memory = "8GB" // 4*2
        time = "5h"

        input:
        file "Aligned.bam" from aligned_bam

        output:
        file "Aligned_Sorted_Filtered.bam" into sam_results

        module "samtools/1.9"

        shell:
        '''
        samtools sort Aligned.bam | samtools view -b -F 1804 -q 2 -o Aligned_Sorted_Filtered.bam
        samtools index Aligned_Sorted_Filtered.bam
        '''

}

process picard_readgroups{

        executor = 'lsf'
        queue = 'premium'
        clusterOptions = '-P acc_apollo -sp 50'
        cpus = 1
        memory = "18GB" // 6*3
        time = "18h"

        input:
	file "Aligned_Sorted_Filtered.bam" from sam_results
        ref_genome

        output:
        file "Aligned_sorted_filtered_RG.bam" into Read_Groups
	set "Aligned_SF_RG_Marked_Dups.bam", "Marked_Dup_Metrics.txt" into Marked_dups

        module "picard/2.18.4:samtools/1.9"

        shell:
        """
        java -jar \$PICARD AddOrReplaceReadGroups I=Aligned_Sorted_Filtered.bam O=Aligned_sorted_filtered_RG.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
        java -jar \$PICARD MarkDuplicates I=Aligned_sorted_filtered_RG.bam O=Aligned_SF_RG_Marked_Dups.bam M=Marked_Dup_Metrics.txt
        samtools index Aligned_SF_RG_Marked_Dups.bam
        """
}

process Variant_Calling{

        executor = 'lsf'
        queue = 'premium'
        clusterOptions = '-P acc_apollo -sp 50'
        cpus = 1
        memory = "20GB" // 10*2
        time = "10h"

        input:
        set "Aligned_SF_RG_Marked_Dups.bam", "Marked_Dup_Metrics.txt" from Marked_dups
        ref_genome

        output:
        file "Raw_Variants.vcf" into Variants

        module "gatk/4.0.1.2"

        shell:
        """
        gatk HaplotypeCaller -I Aligned_SF_RG_Marked_Dups.bam -ERC GVCF -O Raw_Variants.vcf -R ${ref_genome}
        """
}
