FROM debian:latest
MAINTAINER NajlaBioinfo <najlabioinfo.github.io>


ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update && apt-get -qq -y install curl bzip2 vim git sudo less tree libboost-all-dev libstdc++6
RUN apt-get update && apt install -y procps g++ make && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 
RUN apt-get upgrade
RUN apt-get update

ENV PYTHON_VERSION 3.6

# conda and test dependencies
RUN apt-get -qq update && apt-get -qq -y install curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 conda=4.6 \
    && conda install -y python=3 \
    && conda update conda \
    && conda install -y -c anaconda boost \
    && conda install -y -c bioconda nextflow \
    && conda install -y -c bioconda cutadapt \
    && conda install -y -c bioconda fastqc \
    && conda install -y -c bioconda star \
    && conda install -y -c bioconda augustus \
    && conda install -y multiqc -c bioconda -c conda-forge \
    && conda install -y -c bioconda/label/cf201901 salmon \      
           
    && apt-get -qq -y remove curl bzip2 \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
ENV PATH /opt/conda/bin:$PATH
