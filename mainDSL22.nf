#!/usr/bin/env nextflow


/*
********
MODULES
********
*/

/*Quanting*/
include { index; quant } from './quanting.nf'

VERSION="0.1.0"

log.info "======================================================================="
log.info " Nextflow RNA-Seq Pipeline _template   (v${VERSION})  by NajlaBioinfo       "
log.info "======================================================================="


params.help = ""
params.threads = 1
params.output = "./output"
params.data = "./data"
params.username = "someusername"

/*Hardware params*/
threads = params.threads


/* 
 * enables modules 
 */
nextflow.enable.dsl = 2


if(params.help) {
    log.info ''
    log.info 'Check Quality Template'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run . -profile template [options]'
    log.info ''
    log.info 'Script Options: '
    log.info '    --username		STR	user name of the workflow'
    log.info '    --threads		INT	Number of threads to use'
    log.info '    --output		DIR	Directory to write output files'
    log.info '    --data		DIR	Directory to read data files'
    log.info '    --help		BOOL	Display help message'
    log.info ''

    return
}
process greeting {
    echo true

    input:
    params.username
    params.output
    params.data

    output:
    println "Hello: $params.username"
    println "ResultPath: $params.output"
    println "DataPath: $params.data" 
    
    script:
    """
    python $baseDir/bin/argcheck.py ${params.username} $baseDir/${params.data}/fastq $baseDir/${params.output}
    """
}

process sayingbye {
    echo true

    input:
    params.username

    output:
    println "Thank you - $params.username - for using  Nextflow template \nMore infos : bhndevtools@gmail.com"

    script:
    """
    echo "Version used:"
    echo "==============="
    conda --version
    nextflow -v
    fastqc -v
    echo "cutadapt" 
    cutadapt --version
    multiqc --version
    echo "salmon, version  0.8.1" 
    echo "Goodbye"
    """
}


params.trimreads = "$baseDir/${params.output}/trimming/tr_R{1,2}.fastq"
read_pairs2_ch= Channel.fromFilePairs( params.trimreads, checkExists:true )


/* 
 * define the data analysis workflow 
 */
workflow  {
    // workflow implementation
    main:
      greeting()
      index()
      quant(read_pairs2_ch)
      sayingbye()
}
