FROM wwliao/hpp_base:latest
LABEL maintainer="wen-wei.liao@yale.edu"

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

WORKDIR /data
