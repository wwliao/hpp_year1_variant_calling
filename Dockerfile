FROM ubuntu:20.04
LABEL maintainer="Wen-Wei Liao wen-wei.liao@yale.edu"

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

# minimap2 v2.21
WORKDIR /opt/minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.21/minimap2-2.21.tar.bz2 && \
    tar jxf minimap2-2.21.tar.bz2 && \
    rm minimap2-2.21.tar.bz2 && \
    cd minimap2-2.21 && \
    make
ENV PATH="/opt/minimap2/minimap2-2.21:/opt/minimap2/minimap2-2.21/misc:$PATH"

# k8 v0.2.5
WORKDIR /opt/k8
RUN wget https://github.com/attractivechaos/k8/releases/download/0.2.5/k8-0.2.5.tar.bz2 && \
    tar jxf k8-0.2.5.tar.bz2 && \
    rm k8-0.2.5.tar.bz2 && \
    cp k8-0.2.5/k8-`uname -s` k8-0.2.5/k8
ENV PATH="/opt/k8/k8-0.2.5:$PATH"

# pbsv v2.6.2
WORKDIR /opt/pbsv
RUN wget https://github.com/PacificBiosciences/pbsv/releases/download/v2.6.2/pbsv && \
    chmod +x pbsv
ENV PATH="/opt/pbsv:$PATH"

# sniffles v1.0.12b (commit: 4ff6ecb7cfc489ec79b5fc0ec8a9c345786198ac)
WORKDIR /opt/sniffles
RUN git clone https://github.com/fritzsedlazeck/Sniffles.git && \
    cd Sniffles && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make
ENV PATH="/opt/sniffles/Sniffles/bin/sniffles-core-1.0.12:$PATH"

# SVIM 2.0.0
WORKDIR /opt/svim
RUN wget https://github.com/eldariont/svim/archive/refs/tags/v2.0.0.tar.gz && \
    tar zxf v2.0.0.tar.gz && \
    rm v2.0.0.tar.gz && \
    cd svim-2.0.0 && \
    pip install .

# Iris v1.0.4
RUN conda install -y irissv=1.0.4

WORKDIR /data
