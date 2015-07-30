#!/bin/bash
set -e

for d in pybedtools-test pbt-test-py2 pbt-test-py3 ; do
    docker build -t $d $d
done

