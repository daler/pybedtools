FROM continuumio/miniconda
MAINTAINER Ryan Dale <dalerr@niddk.nih.gov>
RUN conda install -c daler \
    matplotlib \
    nose \
    numpydoc \
    pip \
    pandas \
    pyyaml \
    sphinx \
    pysam
RUN conda install -c daler \
    tabix \
    bedtools
WORKDIR /opt/pybedtools
