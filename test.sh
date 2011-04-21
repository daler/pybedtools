python setup.py build_ext -i && \
PATH=$PATH:../bedtools/bin/ nosetests --with-doctest --doctest-extension=.rst README.rst .
