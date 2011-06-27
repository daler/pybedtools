
.. _pip: http://pypi.python.org/pypi/pip

.. _Python Package Index: http://pypi.python.org/pypi

.. _Cython: http://cython.org/

.. _argparse: http://code.google.com/p/argparse/

.. _Python: http://www.python.org/

.. include:: includeme.rst

.. _installation:

Installation
------------

Requirements
~~~~~~~~~~~~
First, make sure you have the required packages for installing :mod:`pybedtools`:

#. BEDTools_. The version is not important, but later versions will have more
   features so it's a good idea to get the latest.  Follow the instructions at
   https://github.com/arq5x/bedtools to install, and make sure the programs are
   on your path. That is, you should be able to call `intersectBed` from any
   directory
#. Python_ 2.5 or greater (Python 3 support is coming soon)
#. Cython_, version 0.14.1 or greater
#. argparse_ if you are running Python < 2.7 (Python 2.7 comes with
   `argparse` already)

Both argparse and Cython can be installed with pip_::

    pip install cython argparse

or `easy_install`::

    easy_install cython argparse


Installing :mod:`pybedtools`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For the latest stable version,

#. download the :mod:`pybedtools` source from
   http://pypi.python.org/pypi/pybedtools

#. Unzip the tarball

#. In the newly created directory, run::

        python setup.py install

   (you may need root privleges to do so)

.. note::

    Due to poor support of Cython and C++ (which is what :mod:`pybedtools`
    uses), it is not currently possible to "`easy_install`" :mod:`pybedtools`
    all in one shot from the command line.  As a result, users must download
    the source and install using the method above.


For the bleeding-edge development version, use the Git repository at
http://github.com/daler/pybedtools.


Testing your installation
~~~~~~~~~~~~~~~~~~~~~~~~~
Quick test
``````````
A quick functional test is to create a new script with the following
contents::

    import pybedtools
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    print a.intersect(b)

If this script is called `test.py`, then running it with `python test.py`
should print out::

    chr1	155	200	feature2	0	+
    chr1	155	200	feature3	0	-
    chr1	900	901	feature4	0	+

Running the test suite
``````````````````````
For more extensive testing, you can run the full test suite which requires
`nose` and `PyYAML` to be installed (both `easy_install`-able).  The test suite
will re-compile the Cython extensions, run unit tests and doctests.  To run the
test suite, use::

    sh test.sh $VERSION

where `$VERSION` is the version of Python you'd like to run the tests with
(e.g., `2.7`).

If you have `sphinx` installed (e.g., via `easy_install sphinx`), you can run
the doctests in this documentation by going to the `docs` directory and
running::

    make doctest


