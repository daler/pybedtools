VERSION=$1
python${VERSION} setup.py build_ext -i && \
PYTHONPATH=$PYTHONPATH:. PATH=$PATH:../bedtools/bin/ nosetests-${VERSION} --with-doctest --doctest-extension=.pyx pybedtools/cbedtools.pyx .
