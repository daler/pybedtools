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
:ref:`Tutorial`.  Before running the examples, you need to import
:mod:`pybedtools`:


::

    >>> from pybedtools import bedtool, cleanup

After running the examples, clean up any intermediate temporary files
with::

    >>> cleanup()

Example 1: Save a BED file of intersections, with track line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example saves a new BED file of intersections between ``a.bed`` and
``b.bed``, adding a track line to the output::

    >>> a = bedtool('a.bed')
    >>> a.intersect('b.bed').saveas('a-and-b.bed', trackline="track name='a and b' color=128,0,0")

Example 2: Intersections for a 3-way Venn diagram
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example gets values for a 3-way Venn diagram of overlaps.  This
demonstrates operator overloading of :class:`bedtool` objects::

    >>> # set up 3 different bedtools
    >>> a = bedtool('a.bed')
    >>> b = bedtool('b.bed')
    >>> c = bedtool('c.bed')
    
    >>> (a-b-c).count()  # unique to a
    >>> (a+b-c).count()  # in a and b, not c
    >>> (a+b+c).count()  # common to all 
    >>> # ... and so on, for all the combinations.

Example 3: Flanking sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example gets the genomic sequences of the 100 bp on either side of
features.

The :meth:`bedtool.slop()` method automatically downloads the
``chromSizes`` table from UCSC for the dm3 genome, but you can pass your
own file using the standard BEDTools ``slop`` argument of ``g``.  Note that
this example assumes you have a local copy of the entire dm3 genome saved
as ``dm3.fa``.

::
    
    >>> # set up bedtool
    >>> mybed = bedtool('in.bed')

    >>> # add 100 bp of "slop" to either side.  genome='dm3' tells
    >>> # the slop() method to download the dm3 chromSizes table from
    >>> # UCSC.
    >>> extended_by_100 = mybed.slop(genome='dm3', l=100, r=100)

    >>> # Delete the middle of the now-200-bp-bigger features so 
    >>> # all we're left with is the flanking region
    >>> flanking_features = extended_by_100.subtract('in.bed')

    >>> # Assuming you have the dm3 genome on disk as 'dm3.fa', save the
    >>> # sequences as a new file 'flanking.fa'
    >>> seqs = flanking_features.sequence(fi='dm3.fa').save_seqs('flanking.fa')

    >>> # We could have done this all in one line 
    >>> # (this demonstrates "chaining" of bedtool objects)
    >>> bedtool('in.bed').slop(genome='dm3',l=100,r=100).subtract('in.bed').flanking_features.sequence(fi='dm3.fa').save_seqs('flanking.fa')

    

For more, continue on to the :ref:`Tutorial`.

.. _BEDTools: http://github.com/arq5x/bedtools
.. |dl| image:: images/downloads.png
