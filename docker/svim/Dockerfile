FROM wwliao/hpp_minimap2:latest
LABEL maintainer="wen-wei.liao@yale.edu"

# SVIM 2.0.0
WORKDIR /opt/svim
RUN wget https://github.com/eldariont/svim/archive/refs/tags/v2.0.0.tar.gz && \
    tar zxf v2.0.0.tar.gz && \
    rm v2.0.0.tar.gz && \
    cd svim-2.0.0 && \
    pip install .

WORKDIR /data
