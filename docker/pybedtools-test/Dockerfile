FROM ubuntu:14.04
MAINTAINER Ryan Dale <dalerr@niddk.nih.gov>
RUN apt-get update && apt-get install -qqy \
    build-essential \
    bzip2 \
    ca-certificates \
    git \
    libsm6 \
    libglib2.0-0 \
    xvfb \
    libxext6 \
    libxrender1 \
    wget \
    zlib1g-dev

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda-3.10.1-Linux-x86_64.sh && \
    /bin/bash /Miniconda-3.10.1-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda-3.10.1-Linux-x86_64.sh && \
    /opt/conda/bin/conda install --yes conda pip

ENV PATH /opt/conda/bin:$PATH
RUN git clone https://github.com/arq5x/bedtools2.git /opt/bedtools
RUN cd /opt/bedtools && git checkout v2.24.0 && make
ENV PATH /opt/bedtools/bin:$PATH

ENV DISPLAY :99.0
WORKDIR /opt/tabix
RUN wget http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
RUN tar -xjvf tabix-0.2.6.tar.bz2
RUN cd tabix-0.2.6 && make
ENV PATH /opt/tabix/tabix-0.2.6:$PATH
WORKDIR /
RUN conda update -q conda
RUN conda config --set always_yes yes --set changeps1 no
