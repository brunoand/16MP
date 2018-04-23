# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
# Modified by Christian Diener
# Changes distributed under the MIT License.
FROM debian:jessie

MAINTAINER Christian Diener <mail@cdiener.com>

USER root

# Install all OS dependencies for fully functional notebook server
ENV DEBIAN_FRONTEND noninteractive
RUN echo "deb http://http.debian.net/debian jessie-backports main" > /etc/apt/sources.list.d/jessie-backports.list





    
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
    && apt-get clean
RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen


# Configure environment
ENV CONDA_DIR /opt/conda
ENV PATH $CONDA_DIR/bin:$PATH
ENV SHELL /bin/bash
ENV NB_USER bruno
ENV NB_UID 1000
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

# Create jovyan user with UID=1000 and in the 'users' group
RUN useradd -m -s /bin/bash -N -u $NB_UID $NB_USER && \
    mkdir -p /opt/conda && \
    chown bruno /opt/conda


# Setup Bruno home directory
RUN mkdir /home/$NB_USER/work && \
    mkdir /home/$NB_USER/.local

RUN cd /tmp && \
    mkdir -p $CONDA_DIR && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && \
    /bin/bash Miniconda-latest-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm Miniconda-latest-Linux-x86_64.sh && \
    $CONDA_DIR/bin/conda install --yes conda

# Install Jupyter notebook as jovyan
RUN conda install --yes --quiet \
    numpy \
    matplotlib=1.4.3 \
    pandas \
    scipy \
    r\
    && conda clean -yt
#Install R libraries
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("dada2")' && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("ggplot2")' && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("ShortRead")'

RUN apt install -y --force-yes libsm6 libxext6 libxrender-dev
RUN pip install --upgrade pip
RUN pip install IPython==5.0
RUN pip install wheel
RUN pip install qiime

RUN apt-get install -t jessie-backports openjdk-8-jre-headless --yes

RUN curl -fsSL get.nextflow.io | bash && \
        mv nextflow /opt/
RUN sed -ri 's/backend      : Qt4Agg/backend      : agg/g' /opt/conda/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc
RUN Rscript -e "install.packages('optparse', repos = 'http://cran.us.r-project.org')"
RUN pip install scikit-bio
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("DESeq2")'
RUN apt-get install libssl-dev --yes
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("rhdf5")'
RUN Rscript -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')" && \
    Rscript -e "install.packages('RJSONIO', repos = 'http://cran.us.r-project.org')" && \
    Rscript -e "devtools::install_github('biom', 'joey711')"

#Install the last two software, namely FastQC and BBmap
RUN wget -q -O fastqc.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc.zip -d /opt/ && \
    chmod 755 /opt/FastQC/fastqc && \
    rm fastqc.zip

RUN wget -q -O BBMap_37.09.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_37.09.tar.gz && \
	tar -xzf BBMap_37.09.tar.gz && \
	mv bbmap/ /opt/ && \
	rm BBMap_37.09.tar.gz

