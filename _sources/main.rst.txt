.. include:: includeme.rst

.. _installation:

Installation
------------
:mod:`pybedtools` is a Python package that wraps BEDTools, so you'll need both
installed.

.. _condainstall:

Quick install via `conda`
~~~~~~~~~~~~~~~~~~~~~~~~~
If you're usng the `Anaconda Python distribution
<http://continuum.io/downloads>`_ on Linux, then the following will install
:mod:`pybedtools`::

    conda install -c bioconda pybedtools

You can also install Tabix and BEDTools via conda::

    conda install -c bioconda bedtools htslib

Otherwise, read on for installation on other platforms and in other
environments.

Required
++++++++
:Python_: version 2.7 or greater (Python 3 is supported). If you're setting up
          Python for the first time, the `Anaconda Python distribution
          <http://continuum.io/downloads>`_ is highly recommended.

:BEDTools_:
    The version is not important, but later versions will have more features so
    it's a good idea to get the latest.  Follow the instructions at
    https://github.com/arq5x/bedtools2 to install, and make sure the programs
    are on your path. That is, you should be able to call `bedtools` from
    any directory


:A C/C++ compiler:
    * **Windows:** Use Cygwin, http://www.cygwin.com.  It is probably easiest to select
      all of the 'Devel" group items to be installed.  In addition, ensure the
      `zlib` items are selected for installation as well (using the search
      funciton in the Cygwin install program).
    * **OSX:** Install Xcode from http://developer.apple.com/xcode/
    * **Linux:** `gcc`, usually already installed; on Ubuntu, install with `sudo apt-get install
      build-essentials`

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

Use the Anaconda channel `daler`::

    conda install -c daler pybedtools

This example installs :mod:`pybedtools` and BEDTools into an isolated
environment called `myenv` running Python 3::

    conda create -n myenv -c daler pybedtools bedtools python=3


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

.. code-block:: bash

    git clone https://github.com/daler/pybedtools.git
    cd pybedtools
    git pull
    python setup.py develop



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
Testing the "current installation" means testing the installation into the
current environment, whether this is the system-wide Python, a virtualenv, or
a conda environment.  It requires some additional packages to be installed::

    pip install -r dev-requirements.txt

Run unit tests::

    nosetests -v

Run doctests::

    (cd docs && make doctest)

Test within isolated conda environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the `condatest.sh` script in the top-level dir of the repo. This script
creates a new isolated conda environment and runs both unit tests and doctests.

To run tests under Python 2::

    ./condatest.sh 2

To run tests under Python 3::

    ./condatest.sh 3

Test within isolated Docker containers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This assumes that `Docker <https://www.docker.com/>`_ is installed.

The following command will build Docker two containers -- one for Python 2 and
one for Python 3 -- starting with the base Ubuntu 14.04 container. The first
time the containers are built it will take some time, but they are cached so
subsequent tests will run quickly. Within each of these containers, unit tests
and doctests are run::

    (cd docker && ./full-test.sh)


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
