FROM wwliao/hpp_base:latest
LABEL maintainer="wen-wei.liao@yale.edu"

USER root

# sniffles v1.0.12b (commit: 4ff6ecb7cfc489ec79b5fc0ec8a9c345786198ac)
WORKDIR /opt/sniffles
RUN git clone https://github.com/fritzsedlazeck/Sniffles.git && \
    cd Sniffles && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make
ENV PATH="/opt/sniffles/Sniffles/bin/sniffles-core-1.0.12:$PATH"

# Iris v1.0.4
RUN conda install -y irissv=1.0.4

WORKDIR /data
