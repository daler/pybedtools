#!/bin/bash

git clone https://github.com/arq5x/bedtools.git
cd bedtools
make
BIN=$(pwd)/bin
export PATH=$PATH:$BIN
echo $PATH

