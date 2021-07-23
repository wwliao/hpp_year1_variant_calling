FROM ubuntu:20.04
LABEL maintainer="wen-wei.liao@yale.edu"

USER root
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y \
        git wget autoconf build-essential zlib1g-dev libbz2-dev libcurl4-gnutls-dev \
	liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev

# miniconda3 (python v3.9)
ENV PATH="/opt/conda/bin:$PATH"
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    echo '. /opt/conda/etc/profile.d/conda.sh' > /etc/profile.d/conda.sh

# cmake v3.21.0
WORKDIR /opt/cmake_install
RUN mkdir /opt/cmake && \
    wget https://github.com/Kitware/CMake/releases/download/v3.21.0/cmake-3.21.0-linux-x86_64.sh && \
    bash cmake-3.21.0-linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake && \
    rm -r /opt/cmake_install

# htslib v1.13
WORKDIR /opt/htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar jxf htslib-1.13.tar.bz2 && \
    rm htslib-1.13.tar.bz2 && \
    cd htslib-1.13 && \
    ./configure && \
    make && \
    make install

# samtools v1.13
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar jxf samtools-1.13.tar.bz2 && \
    rm samtools-1.13.tar.bz2 && \
    cd samtools-1.13 && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure && \
    make && \
    make install

# bcftools v1.13
WORKDIR /opt/bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 && \
    tar jxf bcftools-1.13.tar.bz2 && \
    rm bcftools-1.13.tar.bz2 && \
    cd bcftools-1.13/ && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure && \
    make && \
    make install

WORKDIR /data
