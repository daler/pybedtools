#!/bin/bash
DIR=$1
git clone https://github.com/arq5x/bedtools2.git ${DIR}/bedtools2
cd ${DIR}/bedtools2
make
