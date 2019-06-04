/*
===============================================================
 Dseek pipeline wrapped in Nextflow
===============================================================
*/

def helpMessage() {
    println """
    ===================================
     ${workflow.manifest.name} v${workflow.manifest.version}
    ===================================

    Usage:

    nextflow run Dseek.nf --run [Run_Folder] --genome [Genome] -profile [Profiles]

    Mandatory arguments:
      --run                         Run directory
      --genome                      Reference genome
                                      Human:GRCh38.Gencode.v30
                                      Mouse:GRCm38.Gencode.vM2
      -profile                      Configuration profile to use. [local/chimera/liteQC]
      --singleEnd or --pairedEnd

    Optional arguments:
      Data:
        --rawPath                     Relative path to sample dir to search for FASTQ files ["Raw/Illumina"]
        --outPath                     Relative path to sample dir to store outputs ["Processed/Dseek"]
      
      Pipeline parameters:
        --fastqc                      Enable FASTQC
        --qc                          Enable QC pipeline
        --target                      Provide a target interval list for picard hsmetric [use bait if not provided]
        --bait                        Provide a bait interval list for picard hsmetric 
        --filter                      Enable filtering bams
        --peakcalling                 Enable peak calling for chipseq samples

      Nextflow:
        -w/-work-dir                  Working dir ["work"]
        -resume                       Resuming previous run

    """.stripIndent()
    if (params.help_string) {
        println params.help_string.stripIndent()
    }
}

// Show help message
if (params.help || !params.run){
    helpMessage()
    exit 0
}

if (params.run.startsWith("/")) {
    run_path = params.run
    run_name = params.run.replaceAll("(.*/)", "")
} else {
    run_path = "$PWD/${params.run}"
    run_name = run_path.replaceAll("(.*/)", "")
}

/*
 * Create a channel for input read files
 */

if (params.rawPath==".") {
    raw_search_path = "${params.run}/**[._]R{1,2}[._]*fastq.gz"
} else {
    raw_search_path = "${params.run}/${params.rawPath}/**[._]R{1,2}[._]*fastq.gz"
}

if (params.singleEnd){
    Channel
        .fromPath(raw_search_path)
        .ifEmpty { exit 1, "${params.run}/${params.rawPath} - no fastq.gz file detected." }
        .toSortedList()
        .flatten()
        .map { it ->
                   def fastq = it
                   while(it.getParent().toString() != "${params.run}/${params.rawPath}") {it = it.getParent()}
                   def sample_fastq = [it.getName().replaceAll("(.fastq.gz)",""), fastq]
             }
        .groupTuple()
        .set { sample_fastq }
}

if (params.pairedEnd){
    Channel
        .fromFilePairs(raw_search_path)
        .ifEmpty { error "Cannot find any reads matching: ${params.rawPath}/**[._]R{1,2}[._]*fastq.gz" }
        .set {sample_fastq}
}

sample_fastq.into{sample_fastq; sample_fastq_count}
sample_count = sample_fastq_count.count().get()
println "Sample Detected: ${sample_count}"


if (params.bait){
    if (params.target){
        params.Target = params.target
    }
    else {
        params.Target = params.bait
    }
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "Genome: ${params.genome} is not available in the configuration file. Available genomes: [${params.genomes.keySet().join(", ")}]"
}
// Reference index path loaded from configuration
params.fasta = params.genome ? params.genomes[ params.genome ].fasta_file ?: false : false

if (!params.fasta ) {
    exit 1, "Incomplete reference files."
}

println "========================================="
println "Run Folder: $run_path"
println "Genome: ${params.fasta}"
println "Bait: ${params.bait}"
println "Target: ${params.Target}"

/*
 * STEP 1 - align with BWA
 */
process alignment {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/Bams", mode: 'symlink'

    input:
    params.fasta   
    set val(sample_name), file(reads) from sample_fastq

    output:
    set val(sample_name), file("*Aligned.out.bam") into unsorted_bam

    script:
    """
    bwa mem ${params.fasta} ${reads} -t ${task.cpus} | samtools view -@ ${task.cpus} -S -b -o ${sample_name}.Aligned.out.bam
    """
}
unsorted_bam.into { fastqc_bam; bam_to_sort }

/*
 * STEP 2 - FASTQC
 */
process fastqc {
    tag "${sample_name}"
    publishDir "${run_path}/${sample_name}/${params.outPath}/FastQC", mode: 'copy'

    input:
    set val(sample_name), file(reads) from fastqc_bam

    output:
    file "*_fastqc.zip"
    file "*_fastqc.html" into fastqc_report

    when:
    params.fastqc

    script:
    """
    fastqc -q $reads
    """
}

/*
 * STEP 3 - Sort Bam
 */
process sort_bam {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/Bams", mode: 'symlink'
    
    input:
    set val(sample_name), file(bam) from bam_to_sort

    output:
    set val(sample_name), file("*.Aligned.SortedByCoord.bam") into sorted_bam
   
    script:
    """
    samtools sort -@ ${task.cpus} -T temp-samtool-sort -O BAM -o ${sample_name}.Aligned.SortedByCoord.bam ${bam}
    """
}

sorted_bam.into {tobe_annotated_bam; tobe_filtered_bam}

process filter_bam {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/Bams", mode: 'symlink'

    input:
    set val(sample_name), file(bam) from tobe_filtered_bam

    output:
    set val(sample_name), file("*.Aligned.Sorted.Filtered.bam") into filtered_bam

    when:
    params.filter

    script:
    """
    samtools view -@ ${task.cpus} -b -F 1804 -q 2 -O BAM -o ${sample_name}.Aligned.Sorted.Filtered.bam ${bam}
    """
}

/*
 * STEP 4 -  Annotate Read Groups
 */
process annotate_readgroups {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/Bams", mode: 'symlink'

    input:
    set val(sample_name), file(bam) from tobe_annotated_bam

    output:
    set val(sample_name), file("*.Sorted.RG.bam") into RG_bam
    file("*.Sorted.RG.*bai") into RG_bai
    
    script:
    """
    java -jar \${PICARD} AddOrReplaceReadGroups I=${bam} O=${sample_name}.Sorted.RG.bam RGID=${sample_name} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${sample_name}
    samtools index -@ ${task.cpus} ${sample_name}.Sorted.RG.bam
    """
}


/*
 * STEP 5 -  MarkDuplicates
 */
process markdup {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/Bams", mode: 'copy'

    input:
    set val(sample_name), file(bam) from RG_bam
    file(bai) from RG_bai

    output:
    set val(sample_name), file("${sample_name}.Sorted.MarkedDups.bam") into qc_bam
    file ("*.DupMetrics.txt") into dup_metrics
    file("${sample_name}.Sorted.MarkedDups.*bai") into qc_bai
    
    script:
    """
    java -jar -Xmx55g \${PICARD} MarkDuplicates I=${bam} O=${sample_name}.Sorted.MarkedDups.bam M=${sample_name}.DupMetrics.txt
    samtools index -@ ${task.cpus} ${sample_name}.Sorted.MarkedDups.bam 
    """
}

/*
 * STEP 6
 */

qc_bam.into {call_peaks_bam; hs_metrics_bam; insert_size_metrics_bam; wgs_metrics_bam; alignment_metrics_bam}
qc_bai.into {call_peaks_bai; hs_metrics_bai; insert_size_metrics_bai; wgs_metrics_bai; alignment_metrics_bai}

process picard_hs_metrics {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/qc_metrics", mode: 'copy'

    input:
    params.fasta
    params.bait
    params.Target
    set val(sample_name), file(bam) from hs_metrics_bam
    file(bai) from hs_metrics_bai

    output:
    set val(sample_name), file("*.hsMetrics") into hs_metrics

    when: 
    params.qc
    params.bait
   
    script:
    """
    java -jar \$PICARD CollectHsMetrics \\
        I=${bam} O=${sample_name}.hsMetrics \\
        R=${params.fasta} \\
        TARGET_INTERVALS=${params.Target} \\
        BAIT_INTERVALS=${params.bait}
    """
}

process picard_alignment_metrics {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/qc_metrics", mode: 'copy'

    input:
    params.fasta
    set val(sample_name), file(bam) from alignment_metrics_bam
    file(bai) from alignment_metrics_bai

    output:
    set val(sample_name), file("*.AlignmentMetrics") into alignment_metrics

    when:
    params.qc
    
    script:
    """
    java -jar \$PICARD CollectAlignmentSummaryMetrics \\
        I=${bam} \\
        O=${sample_name}.AlignmentMetrics \\
        R=${params.fasta}
    """
}

process picard_insert_size_metrics {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/qc_metrics", mode: 'copy'

    input:
    params.fasta
    set val(sample_name), file(bam) from insert_size_metrics_bam
    file(bai) from insert_size_metrics_bai

    output:
    set val(sample_name), file("*.InsertSizeMetrics") into insert_size_metrics
    file("*.InsertSize_histogram.pdf")

    when:
    params.qc
    params.pairedEnd

    script:
    """
    java -jar \$PICARD CollectInsertSizeMetrics \\
        I=${bam} \\
        O=${sample_name}.InsertSizeMetrics \\
        H=${sample_name}.InsertSize_histogram.pdf VALIDATION_STRINGENCY=SILENT
    """
}

process picard_wgs_metrics {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/qc_metrics", mode: 'copy'

    input:
    params.fasta
    set val(sample_name), file(bam) from wgs_metrics_bam
    file(bai) from wgs_metrics_bai

    output:
    set val(sample_name), file("*.wgsMetrics") into wgs_metrics

    when:
    params.qc

    script:
    """
    java -jar \$PICARD CollectWgsMetrics \\
        I=${bam} \\
        O=${sample_name}.wgsMetrics \\
        R=${params.fasta}
    """
}

/*
 * Merging Metrics
 */
fastqc_report.mix(
    wgs_metrics,
    alignment_metrics,
    insert_size_metrics,
    hs_metrics,
    dup_metrics
    ).set{ all_metrics }

process multiqc {
    publishDir "${run_path}/", mode: 'copy'

    input:
    file(dna_files) from all_metrics.collect()

    output:
    file("multiqc_report.html")

    when:
    params.qc

    script:
    """
    multiqc .
    """
}
process peak_calling {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/MACS_peakcalling", mode: 'copy'
    
    input:
    params.fasta
    set val(sample_name), file(bam) from call_peaks_bam
    file(bai) from call_peaks_bai
    
    output:
    set val(sample_name), file("*")
  
    when:
    params.peakcalling
    
    script:
    """
    macs2 callpeak -t ${bam} -f BAM -g hs -n ${sample_name} -B -q 0.01
    """
}

