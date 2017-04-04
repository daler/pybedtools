#!/bin/bash

# Installs pybedtools and requirements into a fresh Python 2 or 3 environment
# and runs tests.
#
# Note that this script needs to be called from an environment with Cython
# since this does a clean/sdist operation which will Cythonize the source

set -e

PY_VERSION=$1

usage="Usage: $0 py_version[2|3]"
: ${PY_VERSION:?$usage}

log () {
    echo
    echo "[`date`] TEST HARNESS: $1"
    echo
}

log "removing existing env pbtpy${PY_VERSION}"
name=pbtpy${PY_VERSION}
conda env list | grep -q $name && conda env remove -y -n $name

log "starting with basic environment"
conda create -y -n $name --channel bioconda python=${PY_VERSION} \
    bedtools \
    "htslib<1.4" \
    ucsc-bedgraphtobigwig \
    ucsc-bigwigtobedgraph
source activate $name

log "temporarily install cython"
conda install cython

log "force re-cythonizing"
rm -rf dist build
python setup.py clean
python setup.py build
python setup.py sdist

log "uninstall cython"
conda remove cython

log "test installation of sdist"
set -x
(cd dist && pip install pybedtools-*.tar.gz && python -c 'import pybedtools')
set +x

python setup.py clean

log "install test requirements"
source deactivate
conda env list | grep -q $name && conda env remove -y -n $name
conda create -y -n $name --channel bioconda python=${PY_VERSION} \
    --file "requirements.txt" \
    --file "test-requirements.txt" \
    --file "optional-requirements.txt"

source activate $name

log "install pybedtools from setup.py in develop mode to trigger re-cythonizing"
python setup.py develop

log "run tests"
nosetests
(cd docs && make clean && make doctest)
