``pybedtools`` overview and examples
====================================

Python wrapper for Aaron Quinlan's ``BEDtools`` (http://code.google.com/p/bedtools/).


Quick examples for the impatient
-------------------------------

Get the sequences of the 100 bp on either side of features (with automatic
download of chromSizes from UCSC for dm3 genome)::

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

General usage
-------------
The methods in the ``bedtool`` class mimic the command line usage of
``BEDtools`` as closely as possible.  Any flag that can be passed to the
``BEDtools`` programs can be passed as a keyword argument (see examples
below).  Typically, for convenience ``-i``, ``-a``, and ``-b`` are already
passed for you although you can override this by passing these keyword
arguments explicitly. ``BEDtools`` flags that require an argument, for
example the window size ``-w`` in ``windowBed``, are simply passed as
integers.  Boolean flags (e.g., ``-sw`` in ``windowBed``) are passed as
Python booleans.

Typical workflow is to set up a bedtool object with a bed file you already have::

    a = bedtool('in.bed')
    
`a` now references the file ``in.bed`` on disk.  

Using BEDtools from the command line, in order to get a new bed file of the
intersection of this file with another bed file, ``other.bed`` only
returning uniquely intersecting features from ``in.bed``  we might use::

    intersectBed -a in.bed -b other.bed -u > intersection.bed

Using ``pybedtools``::

    b = a.intersect('other.bed', u=True)

This creates a new temp file in ``/tmp`` by default, but you can change where
temp files are saved using ``pybedtools.set_tempdir()``.  To save a file
explicitly and add a trackline that will label it in a genome browser, use ``saveas()``::
    
    b.saveas('intersection.bed',trackline='track name="intersection"')


Most methods return new ``bedtool`` objects, allowing you to chain things together.  To give you 
a flavor of this, here's how you would get 10 random centers of features that are unique to the file
``other.bed``::

    b = a.intersect('other.bed',v=True).feature_centers(100).random_subset(10)
 
When you're done, it's a good idea to clean up temporary files.  Temp files
are never deleted until you explicitly say so::

    pybedtools.cleanup()

If you forget to call ``pybedtools.cleanup()``, you can always manually delete
the files from your temp dir (typically ``/tmp``).  They are the files that
follow the pattern ``pybedtools.*.tmp``.


The keyword arguments to methods are passed directly to the BEDtools
programs.  On/off switches (e.g., ``-c``, ``-u``, or ``-v`` for ``intersectBed``) are
called with a boolean kwarg; others (like ``-f`` for ``intersectBed``) are
passed in like a normal kwarg.  For example::

    a = bedtool('in.bed')
    a.intersect('other.bed', v=True, f=0.5)

Importantly, the docstring for a method automatically includes the help text of
the original ``BEDtools`` program, so you can check which kwargs you can use
directly from the interpreter.

Examples
~~~~~~~~

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

