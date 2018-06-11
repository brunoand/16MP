# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
# Modified by Christian Diener
# Changes distributed under the MIT License.
FROM ubuntu:16.04

MAINTAINER Christian Diener <mail@cdiener.com>

USER root
    
RUN apt-get update && apt-get install -yq --no-install-recommends \
    git \
    wget \
    build-essential \
    pkg-config \
    ca-certificates \
    bzip2 \
    unzip \
    gfortran \
    tar \
    sudo \
    locales \
    libatlas3-base \
    libfreetype6-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libpcre3 \
    liblzma-dev \
    curl \
    libssl-dev \
    libsm6 \
    r-base \
    python \
    python-pip \
    python3 \
    python3-pip \
    python-setuptools \
    python-dev \
    python-tk \
    libgl1-mesa-glx \
    libsm6 \
    libxext6 \
    libxrender-dev \
    mafft \
    openjdk-8-jdk \
    libpcre3-dev \
    && apt-get clean

RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen

# Install newest version of R from source

RUN wget https://cran.r-project.org/src/base/R-3/R-3.4.0.tar.gz && \
    tar -zxf R-3.4.0.tar.gz && \ 
    cd R-3.4.0 && \
    ./configure --with-x=no --enable-R-shlib=yes --with-readline=no --with-cairo=yes && \
    make && \
    make install && \
    export PATH=~/software/R/R-3.4.0/bin:$PATH


RUN Rscript -e "install.packages('optparse', repos = 'http://cran.us.r-project.org')" && \
    Rscript -e "install.packages(c('devtools', 'curl'), repos = 'http://cran.us.r-project.org')" && \
    Rscript -e "install.packages('RcppParallel', repos = 'http://cran.us.r-project.org')" && \
    Rscript -e "install.packages('RJSONIO', repos = 'http://cran.us.r-project.org')" && \
    Rscript -e "devtools::install_github('biom', 'joey711')"

#Install Bioconductor libraries

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("ggplot2")' && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("ShortRead")' && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("dada2")' && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("DESeq2")' && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("rhdf5")'

#RUN apt-get install -t jessie-backports openjdk-8-jre-headless --yes

RUN curl -fsSL get.nextflow.io | bash && \
        mv nextflow /opt/
#RUN sed -ri 's/backend      : Qt4Agg/backend      : agg/g' /opt/conda/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc

#Install the last two software, namely FastQC and BBmap
RUN wget -q -O fastqc.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc.zip -d /opt/ && \
    chmod 755 /opt/FastQC/fastqc && \
    rm fastqc.zip

RUN wget -q -O BBMap_37.09.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_37.09.tar.gz && \
	tar -xzf BBMap_37.09.tar.gz && \
	mv bbmap/ /opt/ && \
	rm BBMap_37.09.tar.gz

#Install Qiime
RUN pip install --upgrade pip
RUN pip install Cython
RUN pip install pandas
RUN pip install numpy
RUN pip install scipy
RUN pip install IPython==5.0
RUN pip install wheel
RUN pip install qiime
RUN pip install seaborn

RUN apt install python3-setuptools --yes \
    python3-dev
RUN pip3 install numpy
RUN pip3 install scikit-bio

#Add aditional scripts
RUN git clone https://github.com/brunoand/16MP.git && \
        chmod +x 16MP/Scripts/*.py && \
	echo 'c' && \
        mv 16MP/Scripts/ /opt/
 
RUN apt install python3-tk

ENV PATH $PATH:/opt/::/opt/FastQC/:/opt/bbmap:/opt/Scripts
