
.. _pip: http://pypi.python.org/pypi/pip

.. _Python Package Index: http://pypi.python.org/pypi

.. _Cython: http://cython.org/

.. _argparse: http://code.google.com/p/argparse/

.. include:: includeme.rst

.. _installation:

Installation
------------

Python requirements
~~~~~~~~~~~~~~~~~~~
* Python 2.5 or greater; Python 3 support coming soon
* argparse_ if you are running Python < 2.7 (Python 2.7 comes with
  `argparse` already)
* Cython_ - part of :mod:`pybedtools` is written in Cython_ for speed

Both argparse and Cython can be installed with pip_::

    pip install cython argparse


or `easy_install`::

    easy_install cython argparse


Installing :mod:`pybedtools`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To install the latest version of :mod:`pybedtools` you have 2 options:

**Option 1:** install from source

    * go to http://github.com/daler/pybedtools 
    * click the Downloads link (|dl|)
    * choose either a ``.tar.gz`` or a ``.zip`` file, whatever you're 
      comfortable with
    * unzip into a temporary directory
    * from the command line, run::

            python setup.py install

      (you may need admin rights to do this)

**Option 2:** use pip_ to automatically download the latest stable version
from the `Python Package Index`_::

        sudo pip install --upgrade pybedtools

.. warning:: 

    These docs are written for the latest version on git.  For docs
    specific to the version you have installed, please see the docs
    included with that version

Installing BEDTools_
~~~~~~~~~~~~~~~~~~~~
:mod:`pybedtools` relies heavily on BEDTools_.  To install BEDTools_,

    * follow the instructions at https://github.com/arq5x/bedtools to
      install

    * make sure all its programs are on your path


Testing your installation
~~~~~~~~~~~~~~~~~~~~~~~~~
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

For more extensive testing:

* If you have `nosetest` installed (e.g., via `pip install nose`) you can
  run the test suite with::

    sh test.sh $VERSION

  where `$VERSION` is the version of Python you'd like to run the tests
  with, e.g., `2.7`.

* If you have `sphinx` installed (e.g., via `pip install sphinx`), you can
  run the doctests by going to the `docs` directory and running::

    make doctest

