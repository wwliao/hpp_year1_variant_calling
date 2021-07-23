FROM wwliao/hpp_minimap2:latest
LABEL maintainer="wen-wei.liao@yale.edu"

USER root

# SVIM-asm v1.0.2
WORKDIR /opt/svimasm
RUN wget https://github.com/eldariont/svim-asm/archive/refs/tags/v1.0.2.tar.gz && \
    tar zxf v1.0.2.tar.gz && \
    rm v1.0.2.tar.gz && \
    cd svim-asm-1.0.2 && \
    pip install .

WORKDIR /data
