FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive

# pysam needs gcc and compression libs; everything else is pure Python
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/vatools
COPY setup.py .
COPY vatools/ vatools/

RUN pip install --no-cache-dir .

WORKDIR /data
