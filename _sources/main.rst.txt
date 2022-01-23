.. include:: includeme.rst

.. _installation:

Installation
------------
:mod:`pybedtools` is a Python package that wraps BEDTools, so you'll need both
installed.

.. _condainstall:

Install via `conda`
~~~~~~~~~~~~~~~~~~~
This is by far the easiest option.  If you're usng the `Anaconda Python
distribution <http://continuum.io/downloads>`_ on Linux, then the following
will install :mod:`pybedtools`::

    conda install --channel conda-forge --channel bioconda pybedtools

You can also install Tabix and BEDTools via conda::

    conda install --channel conda-forge --channel bioconda bedtools htslib

Otherwise, read on for installation on other platforms and in other
environments.

Required
++++++++
:Python_: version 3.6 or greater (Python 3 is supported). If you're setting up
          Python for the first time, the `Anaconda Python distribution
          <http://continuum.io/downloads>`_ is highly recommended.

:BEDTools_:
    The version is not important, but later versions will have more features so
    it's a good idea to get the latest.  Follow the instructions at
    https://github.com/arq5x/bedtools2 to install, and make sure the programs
    are on your path. That is, you should be able to call `bedtools` from
    any directory


:A C/C++ compiler:
    * **OSX:** Install Xcode from http://developer.apple.com/xcode/
    * **Linux:** `gcc`, usually already installed; on Ubuntu, install with `sudo apt-get install
      build-essentials`
    * **Windows:** may work with conda compliers or Cygwin but this is
      untested. Windows is not supported.

Optional
++++++++
The following external tools are **optional**:

:Tabix_ [`download page`_]:
    Required for fast, random access to BED/GFF/GTF/VCF files by providing
    a chrom:start-stop coordinate.  Similar to the above, you should be able to
    call `tabix` from any directory.


Installing :mod:`pybedtools`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Install latest release via `conda` (recommended)
++++++++++++++++++++++++++++++++++++++++++++++++

See :ref:`condainstall` section above.


Install latest release using `pip`
++++++++++++++++++++++++++++++++++

:mod:`pybedtools` is on PyPI, so you can install via `pip` like most Python
packages. Depending on your Python installation, this may require admin
rights::

    pip install pybedtools


Install development version via github
++++++++++++++++++++++++++++++++++++++

Assumptions:

1. `git` is installed
2. Cython is installed (`conda install cython` or `pip install cython`)


The following commands will clone the repository
.. code-block:: bash

    git clone https://github.com/daler/pybedtools.git
    cd pybedtools

The only time the C++ files will be rebuilt from Cython .pyx source is if the
`cythonize` subcommand is used. To rebuild the C++ files using Cython, run:

.. code-block:: bash

    python setup.py cythonize

To install in develop mode, where changes to Python files will be picked up
without having to re-install, use:

.. code-block:: bash

    python setup.py develop

The above will not update when the .pyx files are updated, so if the Cython
source files have been changed, run:

.. code-block:: bash

    python setup.py cythonize develop


See `python setup.py --usage` for more information.


Quick test
~~~~~~~~~~
Paste the following into a new file called `mytest.py`::

    import pybedtools
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    print a.intersect(b)

Run the script with `python mytest.py`. You should get the results::

    chr1	155	200	feature2	0	+
    chr1	155	200	feature3	0	-
    chr1	900	901	feature4	0	+


Running tests, compiling docs
-----------------------------

There are several modes of testing described below, and in each mode both unit
tests and doctests can be run.

The following instructions assume that you have a working copy of the
:mod:`pybedtools` repository and that you're in the top-level dir of repo,
e.g., by running::

    git clone https://github.com/daler/pybedtools.git
    cd pybedtools


Test current installation
~~~~~~~~~~~~~~~~~~~~~~~~~
To test within the existing installation, install the additional packages for
testing::

    conda install --channel conda-forge --channel bioconda \
      --file requirements.txt \
      --file test-requirements.txt \
      --file optional-requirements.txt

Then run unit tests along with module doctests::

    pytest --doctest-modules

Finally, run sphinx doctests::

    (cd docs && make doctest)

Test within isolated conda environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the `condatest.sh` script in the top-level dir of the repo. This script
creates a new isolated conda environment and runs both unit tests and doctests.

To run tests under Python 2::

    ./condatest.sh 2

To run tests under Python 3::

    ./condatest.sh 3

Compile docs
~~~~~~~~~~~~
To compile the docs, from the top-level `pybedtools` directory::

    (cd docs && make html)


Then point a browser to `docs/build/html/index.html`.

Contributing
~~~~~~~~~~~~
Any and all contributions are welcome.  Here's how to contribute:

#. Fork the `pybedtools repository <http://github.com/daler/pybedtools>`_ on
   github (see `forking help <http://help.github.com/fork-a-repo/>`_).

#. Make your changes/fixes/improvements locally.

#. Optional, but much-appreciated: write some tests for your changes.
        (Don't worry about integrating your tests into the test framework. You
        can just attach the tests either as a commited script or as comments to
        the commit and I can integrate them later)

#. Send a pull request (see `pull request help <http://help.github.com/send-pull-requests/>`_)
