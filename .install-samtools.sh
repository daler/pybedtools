#!/bin/bash

VERSION=$1
DIR=$2
cd ${DIR}
wget http://sourceforge.net/projects/samtools/files/samtools/${VERSION}/samtools-${VERSION}.tar.bz2
tar -xjvf samtools-${VERSION}.tar.bz2
cd samtools-${VERSION}
make
