#!/bin/bash

set -e

conda env remove -y -n deleteme
conda create -n deleteme python=2
source activate deleteme
conda remove setuptools
python setup.py install
(cd && python -c 'import pybedtools; print pybedtools.__file__')
source deactivate
