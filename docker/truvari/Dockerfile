FROM wwliao/hpp_base:latest
LABEL maintainer="wen-wei.liao@yale.edu"

USER root

# truvari v3.2.0-dev (commit: 0e43cb8c878baeb89f23576f3c86a95420364afb)
WORKDIR /opt/truvari
RUN git clone https://github.com/ACEnglish/truvari.git && \
    cd truvari && \
    conda install python=3.8 && \
    python -m pip install --upgrade pip setuptools wheel && \
    python setup.py sdist bdist_wheel && \
    pip install dist/Truvari-3.2.0.dev0.tar.gz

WORKDIR /data
