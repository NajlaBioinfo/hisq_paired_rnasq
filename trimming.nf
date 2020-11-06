#!/usr/bin/env nextflow


VERSION="0.1.0"

log.info "========================================"
log.info " Nextflow RNA-Seq Trimming Module  "
log.info "========================================"


/* 
 * Cutadapt parameter
 * Illumina HiSeq Paired 
 * RNA-Seq data
*/
params.read1 = "$baseDir/${params.output}/data/fastq/SRR12687479_300_1.fastq"
params.read2 = "$baseDir/${params.output}/data/fastq/SRR12687479_300_2.fastq"

params.trimoutput = "$baseDir/${params.output}/trimming"
params.checkquality="$baseDir/${params.output}/qualitycheck"

params.adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.adapter3="AGATCGGAAGAG"
params.adapter4="CTCTTCCGATCT"
params.read1trimed1="tr1_R1.fastq"
params.read1trimed2="tr2_R1.fastq"
params.read1trimed3="tr_R1.fastq"
params.read2trimed1="tr1_R2.fastq"
params.read2trimed2="tr2_R2.fastq"
params.read2trimed3="tr_R2.fastq"
process cutadapt {
    publishDir params.trimoutput, mode:'copy'

    input:
    params.adapter1
    params.trimoutput
    params.checkquality
    params.read1
    params.read2
    params.read1trimed1
    params.read1trimed2

    output:
    println "Step 2/4 Quality Check = Cutadapt"
    println "Cutadapt .... Launched."


    script:
    """
    cutadapt --version
    cutadapt -a ${params.adapter3} -A ${params.adapter3} -o ${params.trimoutput}/${params.read1trimed1} -p ${params.trimoutput}/${params.read2trimed1} ${params.read1} ${params.read2} > ${params.trimoutput}/"report_cutadapt1.txt"
    cutadapt -a ${params.adapter4} -A ${params.adapter4} -o ${params.trimoutput}/${params.read1trimed2} -p ${params.trimoutput}/${params.read2trimed2} ${params.read1} ${params.read2}  > ${params.trimoutput}/"report_cutadapt2.txt"
    cutadapt -a ${params.adapter4} -A ${params.adapter4} -o ${params.trimoutput}/${params.read1trimed3} -p ${params.trimoutput}/${params.read2trimed3} ${params.trimoutput}/${params.read1trimed1} ${params.trimoutput}/${params.read2trimed1} > ${params.trimoutput}/"report_cutadapt3.txt"
    cd ${params.trimoutput}
    mv *.txt ${params.checkquality}/multiqc/
    """
}