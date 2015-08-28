#!/bin/bash

set -e
set -x

# Use Agg backend for matplotlib, which avoids X server errors
mplrc=$(python -c 'from matplotlib import matplotlib_fname as mf; print(mf())')
mkdir -p ~/.config/matplotlib
cp $mplrc ~/.config/matplotlib
sed -i "s/: Qt4Agg/: Agg/g" ~/.config/matplotlib/matplotlibrc

# The repo should have been exported to the container; make a copy and do
# a completely clean installation on that copy before running tests.
cd ~
cp -r /opt/pybedtools .
cd pybedtools
python setup.py clean
python setup.py develop
nosetests
