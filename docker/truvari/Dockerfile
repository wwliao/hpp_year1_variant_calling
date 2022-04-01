FROM wwliao/hpp_base:latest
LABEL maintainer="wen-wei.liao@yale.edu"

USER root

# truvari v3.2.0 (commit: 37b6c3790fc6a96717542635249143be06619c90)
WORKDIR /opt/truvari
RUN git clone -b release/v3.2.0 https://github.com/ACEnglish/truvari.git && \
    cd truvari && \
    conda install python=3.8 && \
    python -m pip install --upgrade pip setuptools wheel && \
    python setup.py sdist bdist_wheel && \
    pip install dist/Truvari-3.2.0.tar.gz

WORKDIR /data
