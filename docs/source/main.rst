Overview
========

Installation
------------
To use :mod:`pybedtools` you'll need the latest version of the package and
the latest version of BEDTools_.

1. To install the latest version of :mod:`pybedtools`:

    * go to http://github.com/daler/pybedtools 
    * click the Downloads link (|dl|)
    * choose either a ``.tar.gz`` or a ``.zip`` file, whatever you're 
      comfortable with
    * unzip into a temporary directory
    * from the command line, run::
            
            python setup.py install

      (you may need admin rights to do this)

2. To install BEDTools_

    * follow the instructions at https://github.com/arq5x/bedtools to install
    * make sure all its programs are on your path


Three brief examples
--------------------
Here are three examples to show typical usage of :mod:`pybedtools`.  More
info can be found in the docstrings of :mod:`pybedtools` methods and in the
:ref:`Tutorial`

Example 1
~~~~~~~~~
Save a new BED file of intersections between ``a.bed`` and
``b.bed``, adding a track line to the output::

    >>> from pybedtools import bedtool
    >>> a = bedtool('a.bed')
    >>> a.intersect('b.bed').saveas('a-and-b.bed', trackline="track name='a and b' color=128,0,0")

Example 2
~~~~~~~~~
Get values for a 3-way Venn diagram of overlaps.  This
demonstrates operator overloading of :class:`bedtool` objects::

    >>> from pybedtools import bedtool
    >>> # set up 3 different bedtools
    >>> a = bedtool('a.bed')
    >>> b = bedtool('b.bed')
    >>> c = bedtool('c.bed')
    
    >>> (a-b-c).count()  # unique to a
    >>> (a+b-c).count()  # in a and b, not c
    >>> (a+b+c).count()  # common to all 
    >>> # ... and so on, for all the combinations.

Example 3
~~~~~~~~~
Get the sequences of the 100 bp on either side of features.

The :meth:`bedtool.slop()` method automatically downloads the
``chromSizes`` table from UCSC for the dm3 genome, but you can pass your
own file using the standard BEDTools ``slop`` argument of ``g``.  Note that
this example assumes you have a local copy of the entire dm3 genome saved
as ``dm3.fa``.

::

    >>> from pybedtools import bedtool
    >>> bedtool('in.bed').slop(genome='dm3',l=100,r=100).subtract('in.bed')
    >>> flanking_features.sequence(fi='dm3.fa').save_seqs('flanking.fa')

    

For more, continue on to the :ref:`Tutorial`.

.. _BEDTools: http://github.com/arq5x/bedtools
.. |dl| image:: images/downloads.png
