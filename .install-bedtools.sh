#!/bin/bash
DIR=$1
git clone https://github.com/arq5x/bedtools.git ${DIR}/bedtools
cd ${DIR}/bedtools
make
