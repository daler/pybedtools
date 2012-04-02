.. include:: includeme.rst

.. _installation:

Installation
------------

Pre-installation requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First, make sure you have the following requirements installed before
attempting to install :mod:`pybedtools`.

Many of these requirements are used by other genomics and bioinformatics
software, so you may already have them installed:

#. BEDTools_
        The version is not important, but later versions will have more
        features so it's a good idea to get the latest.  Follow the
        instructions at https://github.com/arq5x/bedtools to install, and make
        sure the programs are on your path. That is, you should be able to call
        `intersectBed` from any directory

#. Python_ 2.5 or greater (Python 3 support is coming soon)

#. A C/C++ compiler
    * **On Windows:** Use Cygwin, http://www.cygwin.com.  It is probably easiest to select
      all of the 'Devel" group items to be installed.  In addition, ensure the
      `zlib` items are selected for installation as well (using the search
      funciton in the Cygwin install program).
    * **On OSX:** Install Xcode from http://developer.apple.com/xcode/
    * **On Linux:** `gcc`, usually already installed; on Ubuntu, install with `sudo apt-get install
      build-essentials`

#. **Optional.** samtools_ [`download page`_]
        Required for BAM support.   Like BEDTools, the version is not
        important.  You will get a warning if you try to run :mod:`pybedtools`
        functions that require `samtools`.  The `samtools` programs must be
        available on the path, so you should be able to call `samtools view`
        from any directory.

#. **Optional.** Tabix_ [`download page`_].  
        Required for fast, random access to BED/GFF/GTF/VCF files by providing
        a chrom:start-stop coordinate.  Similar to the above, you should be
        able to call `tabix` from any directory.

#. Python modules: these are the modules that :mod:`pybedtools` uses.  All but
   argparse_ are optional.

    * argparse_:
            installed automatically if using Python <2.7 (it comes with Python
            2.7); used for command line scripts like the Venn diagram scripts

    * **Optional.** nose_
            used for automated testing, not necessary for working with
            :mod:`pybedtools`

    * **Optional.** scipy_
            used for computing statistics for randomization procedures

    * **Optional.** matplotlib_
            used plotting code (:mod:`pybedtools.contrib.plotting` and by the
            `venn_mpl.py` script for making a Venn diagram with annotated
            counts. You can use `venn_gchart.py` to use the Google Charts API
            to make a Venn diagram or the :mod:`pybedtools.contrib.venn_maker`
            if you don't want to install `matplotlib`.


#. pip_ or easy_install_:
        these are used for automated installation of Python packages.  If you
        don't already have these installed,  you can use these commands to get
        `pip`::

            $ curl -O http://python-distribute.org/distribute_setup.py
            $ python distribute_setup.py
            $ curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
            $ python get-pip.py

Notes on supported systems
++++++++++++++++++++++++++
:Windows:

    Windows does not support command line programs easily, so **BEDTools and
    pybedtools are only supported via Cygwin on Windows**.

:OSX:

    Should run on any system with the above pre-requisites satisfied.

:Linux:

    Routinely tested on Ubuntu 10.04 and 11.10, but should run on any GNU/Linux
    distribution with the above pre-requisites satisfied.


Installing :mod:`pybedtools`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Simple installation
+++++++++++++++++++
For most users, the latest stable version will be most appropriate.  To
install, use `pip`_ or `easy_install`_ to automatically download the code from
the `Python Package Index`_.  As long as you have pip_ or easy_install_, the
following commands should work on Linux, Windows (under Cygwin_), and OSX::

    pip install pybedtools

or::

    easy_install pybedtools

Done! You can now run a quick test of your installation:

Note that you may need to be root in order to install.  If you do not have root
privileges (e.g., if you are installing in your user directory on a cluster), then
use the following options::

    pip install pybedtools --user

or::

    easy_install pybedtools --prefix /dir_you_can_write_to


Quick test
``````````
A quick functional test is to create a new script with the following
contents, which uses example data shipped with :mod:`pybedtools`::

    import pybedtools
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    print a.intersect(b)

If this script is called `test.py`, then running it with `python test.py`
should print out::

    chr1	155	200	feature2	0	+
    chr1	155	200	feature3	0	-
    chr1	900	901	feature4	0	+


Installation for developers or for running tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Clone the repository::

    git clone git@github.com:daler/pybedtools.git


Developer installation depends on the version.

**Version <= 0.6:**  The first command compiles the Cython code and the second is the typical Python
install command::

    python build.py
    sudo python setup.py develop


**Version > 0.6:** After version 0.6, the developer installation was
streamlined, auto-detecting whether you have Cython installed or not::

    sudo python setup.py develop

The `develop` option means that any changes you make to the code will be
reflected system-wide.  However, if you make changes to any of the `.pyx`
files, you will need to again run the command::


    # make changes to a .pyx file....
    $ sudo python setup.py develop


You should now be able to run the "quick test" above.

More detailed developer installation instructions
+++++++++++++++++++++++++++++++++++++++++++++++++
The very latest code is available from GitHub, at
http://github.com/daler/pybedtools.  You'll need to have `git
<http://git-scm.com/>`_ installed.

#. Get the source by cloning the git repository at
   http://github.com/daler/pybedtools, which will create a new `pybedtools`
   directory

#. Move to the newly created directory

#. *[optional]* Re-compile the extensions by running Cython_.  Typically, only
   developers making changes to the `.pyx`, `.h`, or `.cpp` code will need to
   do this.  Cython needs to be installed for this step.

   ::

        python build.py

#. *[optional]* Run the tests.  nosetests_  and PyYAML_ will be installed if they
   are not already available::

        python setup.py nosetests

        # or to test with a specific version of Python:
        python2.5 setup.py nosetests
        python2.6 setup.py nosetests


#. **Install** :mod:`pybedtools`::

        python setup.py install

#. *[optional]* Install `sphinx`_ if needed (`easy_install sphinx`), then run the
   Sphinx doctests::

        cd docs && make doctest

#. *[optional]* Compile the HTML documentation in the `docs` directory of the
   source tree (this also needs `sphinx`_), then point your browser to
   `docs/build/html/index.html`::

        cd docs && make html

#. *[optional]* Compile the PDF version of the documentation (needs `sphinx_`
   and LaTeX installed), then view `docs/build/latex/pybedtools.pdf` (or the
   copy at `docs/build/html/pybedtools_manual.pdf`)::

        cd docs && make latexpdf

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
