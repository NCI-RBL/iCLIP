# Custom version of iCount which contains an 
# extra option to ultilize multiple threads 
# in demultiplex and cutadapt sub commands:
# https://github.com/skchronicles/iCount
FROM ubuntu:18.04

LABEL maintainer=kuhnsa@nih.gov

# Create Container filesystem specific 
# working directory and opt directories
USER root
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# This section installs system 
# packages required for your 
# project. If you need extra 
# system packages add them here.
# iCount requires argparse 
# python package. Installs 
# python/3.6, which includes 
# argparse by default.
# A child python dependency 
# called 'Pillow', relies
# on the following external
# apt-get dependencies:
# https://pillow.readthedocs.io/en/latest/installation.html#building-on-linux
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        gfortran \
        git \
        g++ \
        make \
        libtiff5-dev \
        libjpeg8-dev \
        libopenjp2-7-dev \
        libfreetype6-dev \
        liblcms2-dev \
        libwebp-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libxcb1-dev \
        libatlas-base-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libssl-dev \
        libcurl4-gnutls-dev \
        liblzma-dev \
        libbz2-dev \
        pigz \
        python \
        python3 \
        python3-dev \
        python3-pip \
        python3-setuptools \
        python3-tk \
        pandoc \
        tcl8.6-dev \
        tk8.6-dev \
        wget \
        zlib1g-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install htslib/1.10
RUN wget https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2 \
    && tar -vxjf htslib-1.10.tar.bz2 \
    && rm htslib-1.10.tar.bz2 \
    && cd htslib-1.10 \
    && make \
    && cd /opt2

# Install samtools/1.10
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 \
    && tar -vxjf samtools-1.10.tar.bz2 \
    && rm samtools-1.10.tar.bz2 \
    && cd samtools-1.10 \
    && make \
    && cd /opt2

# Install BEDTools
# Requires version>=2.26, where
# merge command reports strand
# Also, requires python (already 
# satified) to build from source 
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz \
    && tar -xvf bedtools-2.27.1.tar.gz \
    && rm bedtools-2.27.1.tar.gz \
    && cd bedtools2/ \
    && make \
    && make install \
    && cd /opt2

# Install STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.6.1a.tar.gz \
    && tar -xvzf 2.6.1a.tar.gz \
    && rm 2.6.1a.tar.gz \
    && cd STAR-2.6.1a/source/ \
    && make STAR \
    && cd /opt2

# Install custom version of iCount
# includes multithreading changes
# demultiplex and cutadapt subcmds.
# The output from the custom version 
# and the orignial version should 
# produce identical output. 
# The changes that have been made
# should only effect runtime, (i.e.
# it should be significantly faster).
RUN pip3 install numpy pandas pysam pybedtools numpydoc matplotlib
RUN git clone https://github.com/skchronicles/iCount.git \
    && cd iCount/ \
    && pip3 install .[docs,test] \
    && cd /opt2

# Add Dockerfile to container FS
# and export environment variables
# and ensure proper permissions
ADD Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2 \
    && ln -sf /usr/bin/python3 /usr/bin/python
ENV PATH="$PATH:/opt2/samtools-1.10:/opt2/htslib-1.10:/opt2/STAR-2.6.1a/source"

# Reset working directory
WORKDIR /data2