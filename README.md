Celegansflow
=========================
>>> This is a simple Workflow dedicated to c.elegans Hiseq paired-end RNA-seq data built with nextflow.
>>>>> This workflow treat 5 steps( fastqc before trimming, cutadapt trimming , salmon quanting, fastqc after trimming and multiqc report.)


Prerequisites
-------------
  - Nextflow
  - Java 1.7+
  - Docker

Quickstart
==========
Install Nextflow
----------------
```
$ curl -fsSL get.nextflow.io | bash
$ ./nextflow
$ mv nextflow /usr/local/bin
```

Clone Github Repository
-----------------------
```
$ git clone https://github.com/Najlabioinfo/hisq_paired_rnasq.git
$ cd celegans_nfl
```

Pull Docker Image
------------------
```
$ docker pull najlabioinfo/hiseqcelegans
$ cd /home/najlabioinfo/hisq_paired_rnasq
```
Create Docker Image 
------------------ 
- New docker policies the docker image may disapear but you can build it from the dockerfile.
```
$ cd Docker
$ docker build -t najlabioinfo/hiseqcelegans .
```

Usage
-----
```
$ nextflow run . --help
```
or 
```
$ nextflow run main.nf --username myname --output project/ --data project/data
```
or
```
$ nextflow run . -profile template --threads 2 --output output
```
or for DSL2 nextflow synthax version
```
$ python rnaseq-qualcheck.py testername  project/data project
```

Pipeline Options
----------------
Option | Description
--------- | -----------
help | `Display help message.`
threads | `Number of threads to use for each process.`
output | `Directory to write output files to.`
data | `Reads path`
username | `The name of who launched the pipeline`


References
===========
* https://github.com/nf-core/smrnaseq
* https://github.com/nextflow-io/nfcamp-tutorial


Contact
=======
bhndevtools@gmail.com

License
=======
MIT