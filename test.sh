#!/bin/bash
VERSION=$1
python${VERSION} setup.py nosetests

# clean up from venn diagram script tests
if [ -f out.png ]
then
    rm out.png
fi
