FROM pybedtools-test
MAINTAINER Ryan Dale <dalerr@niddk.nih.gov>
RUN conda install \
    matplotlib \
    nose \
    numpydoc \
    pip \
    pandas \
    python=3 \
    pyyaml \
    sphinx
RUN pip install pysam
WORKDIR /opt/pybedtools
