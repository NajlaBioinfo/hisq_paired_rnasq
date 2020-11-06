#!/usr/bin/env nextflow


VERSION="0.1.0"

log.info "========================================"
log.info " Nextflow RNA-Seq Quanting Module  "
log.info "========================================"

/*
Salmon Index*/
params.transcriptome="$baseDir/${params.output}/data/genomes/PRJNA13758_WBPS15_transcripts.fa"
transcriptome_file = file(params.transcriptome)
params.index_ch="$baseDir/${params.output}/salmon/"

process index {
    tag "salmon index"
    cpus 1
    disk '2G'
    executor 'local'

    input:
    params.index_ch
    params.transcriptome


    output:
    println "Salmon Indexing .... Launched."
    /*file 'index' into index_ch*/

    script:
    """
    salmon index --threads $task.cpus -t ${params.transcriptome} -i ${params.index_ch}/index
    """
}
/*
Salmon Quant
params.trimreads = "$baseDir/${params.output}/trimming/tr_R{1,2}.fastq"
*/
params.indexedtranscript="$baseDir/${params.output}/salmon/index/"
params.indexpairid="$baseDir/${params.output}/salmon/"
params.checkquality="$baseDir/${params.output}/qualitycheck"
params.mqcoutput = "${params.checkquality}/multiqc"

process quant {
    tag "$pair_id"
    
    cpus 1

    input:
    /*file index from index_ch*/
    params.index_ch
    params.mqcoutput
    tuple val(pair_id), path(read_pairs2_ch)

    output:
    println "Salmon Quant .... Launched."
    /*file(pair_id) into quant_ch*/

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i ${params.indexedtranscript} -1 ${read_pairs2_ch[0]} -2 ${read_pairs2_ch[1]} -o ${params.indexpairid}
    cd ${params.indexpairid}
    cp -r * ${params.mqcoutput}/
    """
}