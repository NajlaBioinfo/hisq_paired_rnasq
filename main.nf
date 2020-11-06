#!/usr/bin/env nextflow


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

/* 
 * Fastqc parameter
 */
params.reads = "$baseDir/${params.output}/data/fastq/SRR12687479_300_{1,2}.fastq"
params.checkquality="$baseDir/${params.output}/qualitycheck"
Channel 
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .set { read_pairs_ch }
process fastqc_bt {
    tag "FASTQC on $sample_id Reads Before trimming"

    input:
    set sample_id, file(reads) from read_pairs_ch

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
 * Cuttadapt parameter
 * Illumina HiSeq Paired 
 * RNA-Seq data
*/
params.read1 = "$baseDir/${params.output}/data/fastq/SRR12687479_300_1.fastq"
params.read2 = "$baseDir/${params.output}/data/fastq/SRR12687479_300_2.fastq"
params.trimoutput = "$baseDir/${params.output}/trimming"
params.adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.adapter3="AGATCGGAAGAG"
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
    cutadapt -a ${params.adapter1} -A ${params.adapter1} -o ${params.trimoutput}/${params.read1trimed1} -p ${params.trimoutput}/${params.read2trimed1} ${params.read1} ${params.read2} > ${params.trimoutput}/"report_cutadapt1.txt"
    cutadapt -a ${params.adapter2} -A ${params.adapter2} -o ${params.trimoutput}/${params.read1trimed2} -p ${params.trimoutput}/${params.read2trimed2} ${params.read1} ${params.read2}  > ${params.trimoutput}/"report_cutadapt2.txt"
    cutadapt -a ${params.adapter2} -A ${params.adapter2} -o ${params.trimoutput}/${params.read1trimed3} -p ${params.trimoutput}/${params.read2trimed3} ${params.trimoutput}/${params.read1trimed1} ${params.trimoutput}/${params.read2trimed1} > ${params.trimoutput}/"report_cutadapt3.txt"
    cd ${params.trimoutput}
    mv *.txt ${params.checkquality}/multiqc/
    """
}

/* 
 * Fastqc parameter for trimmed read
 */
params.mqcoutput = "${params.checkquality}/multiqc"
process fastqc_at {
    tag "FASTQC on Trimmed read"

    input:
    params.checkquality
    params.output
    params.trimoutput
    params.mqcoutput
    file('*') from params.trimoutput.collect()
    output:
    println "Step 3/4 Quality Check = FastQC"
    println "Fastqc .... Launched."

    script:
    """
    fastqc -v
    python $baseDir/bin/fastqat.py ${params.trimoutput} ${params.mqcoutput}
    """  
} 

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

/* 
 * Multiqc parameter
 */
 process multiqc {
    publishDir params.mqcoutput, mode:'copy'
       
    input:
    params.trimoutput
    params.checkquality
    file('*') from params.checkquality.collect()

    output:
    println "Step 4/4 Quality Check = MultiQC"
    println "MultiQC .... Launched."
     
    script:
    """
    multiqc --version
    cd ${params.checkquality}/multiqc
    multiqc . 
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
    echo "Goodbye"
    """
}