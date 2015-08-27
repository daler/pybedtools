#1/bin/bash

# Installs pybedtools and requirements into a fresh Python 2 or 3 environment
# and runs tests.
#
# Note that this script needs to be called from an environment with Cython
# since this does a clean/sdist operation which will Cythonize the source

set -e

PY_VERSION=$1

usage="Usage: $0 py_version[2|3]"
: ${PY_VERSION:?$usage}

# Ensure we're starting the environment from scratch
name=pbtpy${PY_VERSION}
conda env remove -y -n $name

# Force the re-Cythonizing
python setup.py clean
python setup.py sdist

conda create \
    -y \
    -c daler \
    -n $name \
    python=${PY_VERSION} \
    bedtools \
    matplotlib \
    sphinx \
    numpydoc \
    tabix \
    pysam \
    nose \
    six \
    pyyaml \
    pandas

source activate $name
pip install -e .
nosetests
(cd docs && make clean && make doctest)
