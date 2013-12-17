#!/bin/bash

VERSION=0.2.6
wget http://sourceforge.net/projects/samtools/files/tabix/tabix-${VERSION}.tar.bz2
tar -xjvf tabix-${VERSION}.tar.bz2
cd tabix-${VERSION}
make
export PATH=$PATH:$(pwd)
echo $PATH
