FROM wwliao/hpp_base:latest
LABEL maintainer="wen-wei.liao@yale.edu"

USER root

# pbsv v2.6.2
WORKDIR /opt/pbsv
RUN wget https://github.com/PacificBiosciences/pbsv/releases/download/v2.6.2/pbsv && \
    chmod +x pbsv
ENV PATH="/opt/pbsv:$PATH"

WORKDIR /data
