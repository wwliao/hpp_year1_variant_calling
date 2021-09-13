FROM wwliao/hpp_base:latest
LABEL maintainer="wen-wei.liao@yale.edu"

USER root

# truvari v3.0.0-dev (commit: badfa6f9a0d18d34f421b2437e05effe8141716a)
WORKDIR /opt/truvari
RUN git clone https://github.com/spiralgenetics/truvari.git && \
    cd truvari && \
    python -m pip install --upgrade pip setuptools wheel && \
    python setup.py sdist bdist_wheel && \
    pip install dist/Truvari-3.0.0.dev0.tar.gz

WORKDIR /data
