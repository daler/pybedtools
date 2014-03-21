.. include:: includeme.rst

.. _installation:

Installation
------------
Notes on supported systems
~~~~~~~~~~~~~~~~~~~~~~~~~~
:Windows:

    Windows does not support command line programs easily, so **BEDTools and
    pybedtools are only supported via Cygwin on Windows**.

:OSX:

    Should run on any system with the below pre-requisites satisfied.

:Linux:

    Routinely tested on Ubuntu 12.04, but should run on any GNU/Linux
    distribution with the below pre-requisites satisfied.

Pre-installation requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Many of these requirements are used by other genomics and bioinformatics
software, so you may already have them installed:

Required
++++++++

:BEDTools_:
    The version is not important, but later versions will have more features so
    it's a good idea to get the latest.  Follow the instructions at
    https://github.com/arq5x/bedtools to install, and make sure the programs
    are on your path. That is, you should be able to call `intersectBed` from
    any directory


:Python_: version 2.5 or greater (Python 3 support is coming soon)

:A C/C++ compiler:
    * **On Windows:** Use Cygwin, http://www.cygwin.com.  It is probably easiest to select
      all of the 'Devel" group items to be installed.  In addition, ensure the
      `zlib` items are selected for installation as well (using the search
      funciton in the Cygwin install program).
    * **On OSX:** Install Xcode from http://developer.apple.com/xcode/
    * **On Linux:** `gcc`, usually already installed; on Ubuntu, install with `sudo apt-get install
      build-essentials`

:Cython_:
    Cython is used to compile C++ and `.pyx` files for :mod:`pybedtools`.


Optional
++++++++
The following external tools are **optional**:


:samtools_ [`download page`_]:
    Required for BAM support.   Like BEDTools, the version is not important.
    You will get a warning if you try to run :mod:`pybedtools` functions that
    require `samtools`.  The `samtools` programs must be available on the path,
    so you should be able to call `samtools view` from any directory.

:Tabix_ [`download page`_]:
    Required for fast, random access to BED/GFF/GTF/VCF files by providing
    a chrom:start-stop coordinate.  Similar to the above, you should be able to
    call `tabix` from any directory.

The following extra Python modules are **optional**:

:nose_:
    used for automated testing, not necessary for working with
    :mod:`pybedtools`

:matplotlib_:
    Used by plotting code (:mod:`pybedtools.contrib.plotting` and by the
    `venn_mpl.py` script for making a Venn diagram with annotated counts. You
    can use `venn_gchart.py` to use the Google Charts API to make a Venn
    diagram or the :mod:`pybedtools.contrib.venn_maker` if you don't want to
    install `matplotlib`.




Installing :mod:`pybedtools`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Root access
+++++++++++
If you have root access to your machine, it's straightforward to install
:mod:`pybedtools`.

From PyPI, with `easy_install`::

    $ sudo easy_install cython
    $ sudo easy_install pybedtools

From PyPI, with `pip`::

    $ sudo pip install cython
    $ sudo pip install pybedtools


Latest development version from git (assumes Cython is installed via one of the
above methods)::

    $ git clone https://github.com/daler/pybedtools.git
    $ cd pybedtools
    $ git pull
    $ sudo python setup.py develop

Then, after any updates on github::

    $ cd pybedtools
    $ git pull
    $ sudo python setup.py develop

No root access
++++++++++++++

If you do not have root access to your machine, there are two alternatives,
Anaconda and virtualenv.

:Anaconda:

    If you do not already have a scientific Python installation, perhaps the
    easiest way to get one is with the `Anaconda distribution
    <http://continuum.io/downloads>`_.  This includes NumPy, matplotlib, and
    all :mod:`pybedtools` dependencies.  Upon installing Anaconda, you can
    install :mod:`pybedtools` with::

        $ pip install pybedtools

:virtualenv:

    The second alternative is to create a virtual environment using the
    `virtualenv` package. The following commands will download, unpack, and
    create a new virtual environment called `myEnv` in your home directory (see
    `virtualenv <https://pypi.python.org/pypi/virtualenv>`_ for more)::

        $ wget https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.10.1.tar.gz
        $ tar -xzvf virtualenv-1.10.1.tar.gz
        $ cd virtualenv-1.10.1/
        $ python virtualenv.py ~/myEnv

    From PyPI, install :mod:`pybedtools` using the new virtualenv `~/myEnv`
    (this also installs Cython)::

        $ ~/myEnv/bin/pip install cython
        $ ~/myEnv/bin/pip install pybedtools


    Latest development version rom git, using the new virtualenv `~/myEnv`
    (this also installs Cython)::

        $ git clone https://github.com/daler/pybedtools.git
        $ cd pybedtools
        $ ~/myEnv/bin/pip install cython
        $ ~/myEnv/bin/python setup.py develop


    From now on, in order to use :mod:`pybedtools`, any time you would normally
    call `python`, instead call `~/myEnv/bin/python`.

Quick test
~~~~~~~~~~
Paste the following into a new file called `mytest.py`::

    import pybedtools
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    print a.intersect(b)

Run the script with `python mytest.py` or, if you used the virtualenv method of
installation described above, run the script with `~/myEnv/bin/python
mytest.py`.  You should get the results::

    chr1	155	200	feature2	0	+
    chr1	155	200	feature3	0	-
    chr1	900	901	feature4	0	+


Installation for developers or for running tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For development work, it's best to use the virtualenv and git method described
above.  The latest code is always available at
https://github.com/daler/pybedtools

The `develop` option means that any changes you make to the code will be
reflected system-wide.  However, if you make changes to any of the `.pyx`
files, you will need to again run the command::

    $ sudo python setup.py develop

This will recompile the Cython files.



.. note::

    :mod:`pybedtools` tests are written using the latest git version of
    BEDTools. You can get the latest git version of BEDTools from
    https://github.com/arq5x/bedtools


Running tests, compiling docs
+++++++++++++++++++++++++++++
There are two sets of tests: unit tests and doctests.  The unit tests need
nosetests_ and PyYAML_ installed, and the doctests need `sphinx`_ and
`numpydoc` installed::

    $ pip install nose
    $ pip install PyYAML
    $ pip install sphinx
    $ pip install numpydoc

For the full tests, you'll also need NumPy and matplotlib (which depends on
NumPy) installed.  These are only used for a small part of :mod:`pybedtools`,
and can be difficult to install if not using the Anaconda distribution
described above.  If you don't have matplotlib, you can ignore the
tests that fail because these libraries can't be found.

.. warning::

    The program `bedtools shuffle` may return different results depending on
    your platform and/or version of C++ stdlib installed (even when using the
    same seed).  The tests for :mod:`pybedtools` are written using Ubuntu 12.04
    with gcc 4.6.3.  Using another platform or stdlib version may result in
    test failures just for tests that use `bedtools shuffle`, similar to those
    reported in `issue #93 <https://github.com/daler/pybedtools/issues/93>`_.

Unit tests
``````````
For the unit tests, in the top-level `pybedtools` directory, run::

    $ sh test.sh

Doctests
````````
For the doctests::

    $ cd docs && make doctests

Compile docs
````````````
To compile the docs, from the top-level `pybedtools` directory::

    $ cd docs && make html


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
