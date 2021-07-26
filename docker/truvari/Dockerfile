FROM wwliao/hpp_base:latest
LABEL maintainer="wen-wei.liao@yale.edu"

USER root

# truvari v2.1.1 (commit: 500fa99546366b07cf849a51a12ce3c348fae735)
WORKDIR /opt/truvari
RUN git clone https://github.com/spiralgenetics/truvari.git && \
    cd truvari && \
    python -m pip install --upgrade pip setuptools wheel && \
    python setup.py sdist bdist_wheel && \
    pip install dist/Truvari-2.1.1.tar.gz

WORKDIR /data
