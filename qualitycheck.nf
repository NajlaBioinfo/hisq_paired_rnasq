#!/usr/bin/env nextflow


VERSION="0.1.0"

log.info "========================================"
log.info " Nextflow RNA-Seq Quality Check Module  "
log.info "========================================"

/* 
 * Fastqc parameter
 
params.reads = "$baseDir/${params.output}/data/fastq/SRR12687479_300_{1,2}.fastq"
read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )*/
params.checkquality="$baseDir/${params.output}/qualitycheck"

process fastqc_bt {
    tag "FASTQC on $sample_id Reads Before trimming"

    input:
    tuple val(sample_id), path(reads)

    output:
    println "Step 1/4 Quality Check = FastQC"
    println "Fastqc .... Launched."

    script:
    """
    fastqc -v
    cd $baseDir/${params.output}/data/fastq/
    fastqc  ${reads}
    mv *fastqc* ${params.checkquality}/multiqc/
    """  
}

/* 
 * Fastqc parameter for trimmed read
 */
params.mqcoutput = "${params.checkquality}/multiqc"
params.trimoutput = "$baseDir/${params.output}/trimming"
process fastqc_at {
    publishDir params.trimoutput, mode:'copy'
    tag "FASTQC on Trimmed read"

    input:
    params.trimoutput
    params.checkquality
    
    output:
    println "Step 1/4 Quality Check = FastQC"
    println "Fastqc .... Launched."

    script:
    """
    fastqc -v
    cd ${params.trimoutput}
    fastqc -q *.fastq
    mv *fastqc* ${params.checkquality}/multiqc/
    """
} 

/* 
 * Multiqc parameter
 */
 process multiqc {
    publishDir params.mqcoutput, mode:'copy'
       
    input:
    params.checkquality


    output:
    println "Step 4/4 Quality Check = MultiQC"
    println "MultiQC .... Launched."
     
    script:
    """
    multiqc --version
    cd ${params.mqcoutput}
    pwd
    python $baseDir/bin/multiqcst.py ${params.mqcoutput}/ ${params.mqcoutput}/ 19
    """
} 