.. include:: includeme.rst

.. _installation:

Installation
------------

Requirements
~~~~~~~~~~~~
First, make sure you have the following requirements installed before
attempting to install :mod:`pybedtools`:

#. BEDTools_. The version is not important, but later versions will have more
   features so it's a good idea to get the latest.  Follow the instructions at
   https://github.com/arq5x/bedtools to install, and make sure the programs are
   on your path. That is, you should be able to call `intersectBed` from any
   directory
#. Python_ 2.5 or greater (Python 3 support is coming soon)
#. A C/C++ compiler
    * Windows: http://www.mingw.org/
    * OSX: Install Xcode from http://developer.apple.com/xcode/
    * Linux: usually already installed; on Ubuntu, install with `sudo apt-get install
      build-essentials`

Installing :mod:`pybedtools`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Simple installation
+++++++++++++++++++
For most users, the latest stable version will be most appropriate.  To
install, use `pip`_ or `easy_install`_ to automatically download the code from
the `Python Package Index`_::

    pip pybedtools

or::

    easy_install pybedtools

You may need to be root in order to install.  If you do not have root
privleges (for example, installing in your user directory on a cluster), then
use the `--prefix` argument to `easy_install` to specify a location where you
have write permission::

    easy_install pybedtools --prefix /dir_you_can_write_to

Done! You can now run a quick test of your installation:

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


Installation for developers or for running tests
++++++++++++++++++++++++++++++++++++++++++++++++
A more flexible way to install :mod:`pybedtools` is as follows:

#. Get the source.  There are two ways to do this:

    #. **Stable version:** download and unzip the :mod:`pybedtools` source from
       http://pypi.python.org/pypi/pybedtools

    #. **Development version:** clone the git repository at
       http://github.com/daler/pybedtools

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

        cd docs && make doctests

#. *[optional]* Compile the HTML documentation (also needs `sphinx`_), then point your browser to
   `docs/build/html.index.html`::

        cd docs && make html

#. *[optional]* Compile the PDF documentation (needs `sphinx_`), then view `manual.pdf`::

        cd docs && make latexpdf
