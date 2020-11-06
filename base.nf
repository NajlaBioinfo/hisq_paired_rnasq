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


process Greeting {
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
    python $baseDir/bin/argcheck.py ${params.username} $baseDir/${params.data} $baseDir/${params.output}
    """
}

params.prefinput = "$baseDir/${params.data}sra/sraIds.txt"
params.fastqinput="$baseDir/${params.data}fastq/"
prefetchIds = Channel.create()
Channel
    .from(file(params.prefinput))
    .splitCsv(sep:'\t', header: true)
    .map { it.Run_s }
    .into(prefetchIds)

process FastqDumpDocker {
    echo true
	publishDir "${params.output}"
	container 'najlabioinfo/celegans:dev'

   input:
    params.prefinput
    params.output
    params.fastqinput
    val idsra from prefetchIds

    script:
    """
    echo "Step 1/3 Fastq-dump ...."
    docker run -t --name $idsra"f"  najlabioinfo/celegans  fastq-dump -I --split-files ${idsra} > $params.fastqinput/
    docker stop $idsra"f"
    docker rm $idsra"f"
    """
}

process SayingBye {
    echo true

    input:
    params.username

    output:
    println "Thank you - $params.username - for using  Nextflow template \nMore infos : bhndevtools@gmail.com"

    script:
    """
    echo "Goodbye"
    """
}

params.prefinput = "$baseDir/${params.data}sra/sraIds.txt"
params.fastqinput="$baseDir/${params.data}fastq/"
prefetchIds = Channel.create()
Channel
    .from(file(params.prefinput))
    .splitCsv(sep:'\t', header: true)
    .map { it.Run_s }
    .into(prefetchIds)

process FastqDumpDocker2 {
    echo true
	publishDir "${params.output}"
	container 'najlabioinfo/condasras:latest'

   input:
    params.prefinput
    params.fastqinput
    val idsra from prefetchIds
    params.output


    script:
    """
    echo "Step 1/3 Fastq-dump ...."
    echo "Version:"
    docker run -t --name "fcontainer"  najlabioinfo/condasras:latest fastq-dump --version

    echo "Stop container ...."
    docker stop fcontainer
    echo "Remove container ...."
    docker rm fcontainer
    echo "Downloading fastq files ...."
    docker run -t --name $idsra"f"   -v /root/ncbi/public/sra -v 271020/data/fastq najlabioinfo/condasras:latest  fastq-dump -X 100 -Z ${idsra}  > $params.fastqinput/$idsra".fastq" 
    echo "Stop the container ...."
    docker stop $idsra"f"
    echo "Remove the container ...."
    docker rm $idsra"f"
    """
}