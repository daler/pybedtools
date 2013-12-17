#!/bin/bash

VERSION=$1
DIR=$2
cd ${DIR}
wget http://sourceforge.net/projects/samtools/files/tabix/tabix-${VERSION}.tar.bz2
tar -xjvf tabix-${VERSION}.tar.bz2
cd tabix-${VERSION}
make
