``pybedtools`` overview and examples
====================================

``pybedtools`` is a Python wrapper for Aaron Quinlan's ``BEDtools`` (http://code.google.com/p/bedtools/).

See full online documentation at http://daler.github.com/pybedtools (this
README contains just installation instructions and a few examples).


Installation
------------

Use ``pip`` to install the latest from the Github repo::

    pip install pybedtools

or ``easy_install``::

    easy_install pybedtools


Quick examples for the impatient
--------------------------------

Get the sequences of the 100 bp on either side of features (with automatic
download of chromSizes from UCSC for dm3 genome).  This assumes you have a
local copy of the entire dm3 genome as ``dm3.fa``::

    import pybedtools
    pybedtools.bedtool('in.bed').slop(genome='dm3',l=100,r=100).subtract('in.bed')
    flanking_features.sequence(fi='dm3.fa').save_seqs('flanking.fa')

Or, get values for a 3-way Venn diagram of overlaps::

    import pybedtools
    a = pybedtools.bedtool('a.bed')
    b = pybedtools.bedtool('b.bed')
    c = pybedtools.bedtool('c.bed')

    (a-b-c).count()  # unique to a
    (a+b-c).count()  # in a and b, not c
    (a+b+c).count()  # common to all 
    # ... and so on, for all the combinations.
    
Intersections, plus adding track names::


    import pybedtools
    a = pybedtools.bedtool('a.bed')
    a.intersect('b.bed').saveas('a-and-b.bed', trackline="track name='a and b'")
   
Creating a :class:`pybedtools.bedtool`:

.. doctest::

    >>> import pybedtools
    >>> a = pybedtools.bedtool('chrX 1 100', from_string=True).saveas('a.bed')

.. doctest::
   :options: +NORMALIZE_WHITESPACE

   >>> b = pybedtools.bedtool('a.bed')
   >>> print b
   chrX     1       100
   <BLANKLINE>




Why use ``pybedtools``?
-----------------------
Using this module allows you to use BEDtools directly from your Python code
without awkward system calls.  The provided ``pybedtools.bedtool`` class
wraps the BEDtools command line programs in an intuitive and easy-to-use
interface.  As a quick illustration of the streamlining possible, here's
how to get the number of features shared between a.bed and b.bed, those
unique to a.bed, and those unique to b.bed::

    from pybedtools import bedtool
    a = bedtool('a.bed')
    b = bedtool('b.bed')
    (a+b).count()    # shared in a and b
    (a-b).count()    # unique to a
    (b-a).count()    # unique to b

In contrast, here's how you'd do the same from the command line:: 

    intersectBed -a a.bed -b b.bed -u | wc -l   # shared in a and b
    intersectBed -a a.bed -b b.bed -u | wc -l   # unique to a
    intersectBed -a b.bed -b a.bed -u | wc -l   # unique to b

To do the same thing in Python, *each* of these lines would have to be
wrapped in awkward, piped ``subprocess.Popen`` calls::
    
    p1 = subprocess.Popen(['intersectBed','-a','a.bed','-b','b.bed','-u'], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['wc','-l'], stdin=subprocess.PIPE)
    results = p2.communicate()[0]
    count = results.split()[-1]

See, ``a+b`` is much easier!

Translation between BEDTools programs and ``pybedtools.bedtool`` methods
------------------------------------------------------------------------

.. currentmodule:: pybedtools

This table summarizes the translation between BEDTools program names
:class:`bedtool` method names.

================= ====================================
BEDTools program  ``pybedtools.bedtool`` method
================= ====================================
`intersectBed`    :meth:`bedtool.intersect`
`subtractBed`     :meth:`bedtool.subtract`
`fastaFromBed`    :meth:`bedtool.sequence`
`slopBed`         :meth:`bedtool.slop`
`mergeBed`        :meth:`bedtool.merge`
`closestBed`      :meth:`bedtool.closest`
`windowBed`       :meth:`bedtool.window`
`groupBy`         :meth:`bedtool.groupBy`
`shuffleBed`      :meth:`bedtool.shuffle`
`sortBed`         :meth:`bedtool.sort`
================= ====================================

There are also methods unique to ``pybedtools`` that provide additional usefulness:

    * :meth:`bedtool.size_filter`
    * :meth:`bedtool.random_subset`
    * :meth:`bedtool.saveas`
    * :meth:`bedtool.cat`
    * :meth:`bedtool.randomstats`
    * :meth:`bedtool.get_genome`
    * :meth:`bedtool.save_seqs`
    * :meth:`bedtool.count`
    * :meth:`bedtool.features`
    * :meth:`bedtool.lengths`

General usage
-------------
.. doctest::

    >>> from pybedtools import bedtool
    >>> a = '''
    ...         chrX 1   100
    ...         chrX 200 500
    ...         chrY 499 600
    ... '''
    >>> b = '''
    ...         chrX 10  60
    ...         chrY 200 500
    ... '''
    >>> a = bedtool(a, from_string=True).saveas('a.bed')
    >>> b = bedtool(b, from_string=True).saveas('b.bed')

Arguments are the same as BEDtools command line programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The methods in the ``bedtool`` class mimic the command line usage of
``BEDtools`` as closely as possible.  Any flag that can be passed to the
``BEDtools`` programs can be passed as a keyword argument in the
corresponding ``pybedtools.bedtool`` method.

On/off switches (e.g., ``-c``, ``-u``, or ``-v`` for ``intersectBed``) are
called with a boolean kwarg; others (like ``-f`` for ``intersectBed``) are
passed in like a normal kwarg with appropriate values.  For example::

    a = bedtool('in.bed')
    a.intersect('other.bed', v=True, f=0.5)

Typically, for convenience ``-i``, ``-a``, and ``-b`` are already
passed for you although you can override this by passing these keyword
arguments explicitly.  The second line above could have equivalently been
called as::

    a.intersect(a='in.bed', b='other.bed', v=True, f=0.5)

Conveniently, the docstring for a method automatically includes the help
text of the original ``BEDtools`` program, so you can check which kwargs
you can use directly from the interpreter.

Typical workflow includes temporary files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Typical workflow is to set up a bedtool object with a bed file you already have::

    a = bedtool('in.bed')
    
`a` now references the file ``in.bed`` on disk.  

Using BEDtools from the command line, we might use the following in order
to get a new bed file of the intersection of this file with another bed
file, ``other.bed`` but only returning uniquely intersecting features from
``in.bed``::

    intersectBed -a in.bed -b other.bed -u > intersection.bed

Using ``pybedtools``::

    b = a.intersect('other.bed', u=True)

This creates a new temp file in ``/tmp`` by default, but you can change
where temp files are saved using ``pybedtools.set_tempdir()``.  To save a
file explicitly and optionally add a trackline that will label it in a
genome browser, use ``saveas()``::
    
    b.saveas('intersection.bed',trackline='track name="intersection"')

Cleaning up temporary files
~~~~~~~~~~~~~~~~~~~~~~~~~~~
When you're done, it's a good idea to clean up temporary files.  Temp files
are never deleted until you explicitly say so::

    pybedtools.cleanup()

If you forget to call ``pybedtools.cleanup()``, you can always manually
delete the files from your temp dir (typically ``/tmp``), though you can
specify this with ``pybedtools.set_tempdir()``.  They are the files that
follow the pattern ``pybedtools.*.tmp``.

Chaining together commands
~~~~~~~~~~~~~~~~~~~~~~~~~~
Most methods return new ``bedtool`` objects, allowing you to chain things
together.  This means that you can chain commands together, just like
piping things together on the command line.  To give you a flavor of this,
here's how you would get 10 random centers of features that are unique to
the file ``other.bed``::

    b = a.intersect('other.bed',v=True).feature_centers(100).random_subset(10)
 

Individual methods
------------------

:meth:`bedtool.intersect`
~~~~~~~~~~~~~~~~~~~~~~~~~

Examples of how to use :meth:`bedtool.intersect`.


First, we set up the files to use:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> # Set up the bedtools
    >>> a = bedtool('a.bed')
    >>> b = bedtool('b.bed')

Here's what they look like:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print a
    chrX	1	100
    chrX	200	500
    chrY    499 600

    >>> print b 
    chrX	10	60
    chrY	200	500

Default is to only report the part that intersects:


.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print a.intersect(b)
    chrX	10	60
    chrY    499 500

Use the `-u` flag of ``intersectBed`` to report full features in :file:`a.bed`
that intersected :file:`b.bed`:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print a.intersect(b, u=True)
    chrX   1    100
    chrY   499  600

Or the opposite -- features in :file:`a.bed` that are not in :file:`b.bed`:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print a.intersect(b, v=True)
    chrX	200	500

Report features in :file:`a.bed`, attaching an additional column indicating how many
features in :file:`b.bed` intersected:

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> print a.intersect(b, c=True)
    chrX	1	100	1
    chrX	200	500	0
    chrY	499	600	1

You can retrieve these counts later using the :meth:`bedtool.counts` method:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> result = a.intersect(b, c=True)
    >>> print result.counts()
    [1, 0, 1]

    

Example: Flanking seqs
----------------------
The ``slop()`` method (which calls ``slopBed``) needs a chromosome size
file.  If you specify a genome name to the ``slop()`` method, it will
retrieve this file for you automatically from the UCSC Genome Browser MySQL
database.

::
    
    import pybedtools
    a = pybedtools.bedtool('in.bed')
    extended = a.slop(genome='dm3',l=100,r=100)
    flanking = extended.subtract(a).saveas('flanking.bed')
    flanking.sequence(fi='dm3.fa')
    flanking.save_seqs('flanking.fa')


Or, as a one-liner::

    pybedtools.bedtool('in.bed').slop(genome='dm3',l=100,r=100).subtract(a).sequence(fi='dm3.fa').save_seqs('flanking.fa')

Don't forget to clean up!::

    pybedtools.cleanup()

Example: Region centers that are fully intergenic
--------------------------------------------------
Useful for, e.g., motif searching::
    
    a = pybedtools.bedtool('in.bed')
    
    # Sort by score
    a = a.sorted(col=5,reverse=True)

    # Exclude some regions
    a = a.subtract('regions-to-exclude.bed')

    # Get 100 bp on either side of center
    a = a.peak_centers(100).saveas('200-bp-peak-centers.bed')


Example: Histogram of feature lengths
-------------------------------------
Note that you need matplotlib installed to plot the histogram.

::

    import pylab as p
    a = pybedtools.bedtool('in.bed')
    p.hist(a.lengths(),bins=50)
    p.show()


Selected examples from http://code.google.com/p/bedtools/wiki/UsageAdvanced
----------------------------------------------------------------------------
Command line::

    intersectBed -a snp.calls.bed -b dbSnp.bed -v | intersectBed -a stdin -b 1KG.bed -v > snp.calls.novel.bed

Same thing, in ``pybedtools``::
    a = pybedtools.bedtool('snp.calls.bed')
    a.intersect('dbSnp.bed',v=True).intersect('1KG.bed',v=True).saveas('snp.calls.novel.bed')
