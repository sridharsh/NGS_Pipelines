# NextGenSeeker

One Paragraph of project description goes here

## Getting Started

NGS data, such as WGS, RNA-Seq, ChIP-Seq, single cell sequencing can generate significant amounts of output data. NextGenSeeker chooses appropriate softwares and tools that are existing for the bioinformatics analysis. The processes involved in the data analysis pipelines for each of the sequencing types includes raw reads quality control, preprocessing, mapping, post-alignment processing, variant calling, followed by variant annotation and prioritization, peakcalling, differential gene expression analysis, read counting and other downstream analysis.

### Quality Controlling the data using FastQC

To assess the quality of raw sequencing data, you will need to load the FastQC tool of version 0.11.8.  
Run FastQC by running the following command:
```
fastqc -f <input_file_format> <read_sample1> <read_sample2>
```

Note: The file format can be either fastq, bam or sam

This step will give a general view on number and length of reads, if there are any contaminating sequences in the sample or low-quality sequences. This step can be done for any of the types of samples before running their respective pipelines.

## D-seek
## R-seek
## ChIP-seek


