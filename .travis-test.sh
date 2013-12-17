#!/bin/bash

set -e

TABIX_VERSION=0.2.6
SAMTOOLS_VERSION=0.1.19
./.install-tabix.sh $TABIX_VERSION ${TRAVIS_BUILD_DIR}
./.install-bedtools.sh ${TRAVIS_BUILD_DIR}
./.install-samtools.sh $SAMTOOLS_VERSION ${TRAVIS_BUILD_DIR}
export PATH=$PATH:${TRAVIS_BUILD_DIR}/tabix-${TABIX_VERSION}
export PATH=$PATH:${TRAVIS_BUILD_DIR}/samtools-${SAMTOOLS_VERSION}
export PATH=$PATH:${TRAVIS_BUILD_DIR}/bedtools/bin
echo $PATH
nosetests -v
cd docs && make clean && make doctest
