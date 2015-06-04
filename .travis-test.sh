#!/bin/bash

set -e

TABIX_VERSION=0.2.6
SAMTOOLS_VERSION=0.1.19
BEDTOOLS_BRANCH=v2.24.0
./.install-tabix.sh $TABIX_VERSION ${TRAVIS_BUILD_DIR}
./.install-bedtools2.sh ${TRAVIS_BUILD_DIR} ${BEDTOOLS_BRANCH}
./.install-samtools.sh $SAMTOOLS_VERSION ${TRAVIS_BUILD_DIR}
export PATH=${TRAVIS_BUILD_DIR}/tabix-${TABIX_VERSION}:$PATH
export PATH=${TRAVIS_BUILD_DIR}/samtools-${SAMTOOLS_VERSION}:$PATH
export PATH=${TRAVIS_BUILD_DIR}/bedtools2/bin:$PATH
echo $PATH
nosetests -v
cd docs && make clean && make doctest

