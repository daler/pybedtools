#!/bin/bash

set -e

echo $PATH

VERSION=$1

python${VERSION} setup.py nosetests -v

# clean up from venn diagram script tests
if [ -f out.png ]
then
    rm out.png
fi
