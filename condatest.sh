#!/bin/bash

set -e

PY_VERSION=$1

usage="Usage: $0 py_version[2|3]"
: ${PY_VERSION:?$usage}


log () {
    echo
    echo "[`date`] TEST HARNESS: $1"
    echo
}



# ----------------------------------------------------------------------------
# sdist and pip install tests
# ----------------------------------------------------------------------------
# Build an environment with just Python and Cython. We do this fresh each time.
log "building fresh environment with just python and cython"
with_cy="pbtpy${PY_VERSION}_sdist_cython"
if conda env list | grep -q $with_cy; then
    conda env remove -y -n $with_cy
fi
conda create -n $with_cy -y --channel conda-forge --channel bioconda python=${PY_VERSION} cython
source activate $with_cy

# Clone the repo -- so we're only catching things committed to git -- into
# a temp dir
log "cloning into temp dir"
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TMP=/tmp/pybedtools-deploy
rm -rf $TMP
git clone $HERE $TMP
cd $TMP

log "cythonizing source files and building source package"
# Cythonize the .pyx filex to .cpp, and build a source package
python setup.py clean cythonize sdist

log "installing source package with pip"
# Install into the environment to verify that everything works (just an import
# test)
(cd dist && pip install pybedtools-*.tar.gz && python -c 'import pybedtools; print(pybedtools.__file__)')


# ----------------------------------------------------------------------------
# Unit tests
# ----------------------------------------------------------------------------
# Deactivate that env, and build another one with all requirements that we'll
# use for unit tests.
source deactivate
no_cy="pbtpy${PY_VERSION}_conda_no_cython"
if ! conda env list | grep -q $no_cy; then
    log "creating environment"

    # pysam not available from bioconda for py37 so remove it from
    # requirements.
    TMPREQS=$(tempfile)
    grep -v pysam requirements.txt > $TMPREQS
    if [[ "$PY_VERSION" == "3.7" ]]; then
        REQS=$TMPREQS
    else
        REQS=requirements.txt
    fi

    # genomepy>=0.8 not available for py27
    TMPOPTREQS=$(tempfile)
    grep -v genomepy optional-requirements.txt > $TMPOPTREQS
    if [[ "$PY_VERSION" == "2.7" ]]; then
        OPTREQS=$TMPOPTREQS
    else
        OPTREQS=optional-requirements.txt
    fi


    conda create -n $no_cy -y \
        --channel conda-forge \
        --channel bioconda \
        python=${PY_VERSION} \
        --file $REQS \
        --file test-requirements.txt \
        --file $OPTREQS
else
    echo "Using existing environment '${no_cy}'"
fi
source activate $no_cy

log "unpacking source package and install with pip install -e into $no_cy env"
mkdir -p /tmp/pybedtools-uncompressed
cd /tmp/pybedtools-uncompressed
tar -xf $TMP/dist/pybedtools-*.tar.gz
cd pybedtools-*
pip install -e .

log "Unit tests"
pytest -v --doctest-modules --ignore pybedtools/test/test_genomepy_integration.py
pytests -v pybedtools/test/test_genomepy_integration.py

# ----------------------------------------------------------------------------
# sphinx doctests
# ----------------------------------------------------------------------------
# Since the docs aren't included in the MANIFEST and therefore aren't included
# in the source distribution, we copy them over from the repo we checked out.
log "copying over docs directory from repo"
cp -r $TMP/docs .

log "sphinx doctests"
(cd docs && make clean doctest)
