python setup.py build_ext -i && \
PATH=$PATH:../bedtools/bin/ nosetests --with-doctest --doctest-extension=.pyx pybedtools/cbedtools.pyx .&& \
# cd docs && make doctest  && make html
