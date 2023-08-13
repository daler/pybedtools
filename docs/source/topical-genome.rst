.. include:: includeme.rst

.. _genomes:

Specifying genomes
==================
This section illustrates the use of genome files for use with BEDTools
programs that need to know chromosome limits to prevent out-of-range
coordinates.

Using BEDTools programs like `slopBed` or `shuffleBed`  from the command
line requires "genome" or "chromsizes" files.  :mod:`pybedtools` comes with
common genome assemblies already set up as a dictionary with chromosomes as
keys and zero-based (start, stop) tuples as values:

.. doctest::

    >>> from pybedtools import genome_registry
    >>> genome_registry.dm3['chr2L']
    (0, 23011544)


The rules for specifying a genome for methods that require a genome are as
follows (use whatever is most convenient):

* Use `g` to specify either a filename or a dictionary
* Use `genome` to specify either an assembly name or a dictionary

Below are examples of each.

As a file
---------

This is the typical way of using BEDTools programs, by specifying an existing genome
file with `g`:

.. doctest::
    :hide:

    >>> fn = pybedtools.chromsizes_to_file(pybedtools.chromsizes('hg19'), fn='hg19.genome')

.. doctest::

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = a.slop(b=100, g='hg19.genome')

.. doctest::
    :hide:

    >>> import os
    >>> os.unlink('hg19.genome')

As a string
-----------
This is probably the most convenient way of specifying a genome.  If the
genome exists in the genome registry it will be used directly.  Alternatively,
if you have `genomepy<https://github.com/vanheeringen-lab/genomepy/>`_ 
installed, you can use the genomepy genome name, such as `hg38`.  In this case, 
the genome file will be located automatically.  Finally, if the genome is not 
in the registry or managed by genomepy, it will automatically be downloaded 
from UCSC.  You must use the `genome` kwarg for this; if you use `g` a string 
will be interpreted as a filename:

.. doctest::

    >>> c = a.slop(b=100, genome='hg19')

As a dictionary
---------------
This is a good way of providing custom coordinates; either `g` or `genome`
will accept a dictionary:

.. doctest::

    >>> d = a.slop(b=100, g={'chr1':(1, 10000)})
    >>> e = a.slop(b=100, genome={'chr1':(1,100000)})

Make sure that all these different methods return the same results

.. doctest::

    >>> b == c == d == e
    True

Converting to a file
--------------------
Since BEDTools programs operate on files, the fastest choice will be to
use an existing file.  While the time to convert a dictionary to a file is
extremely small, over 1000's of files (e.g., for Monte Carlo simulations),
the time may add up.  The function :func:`pybedtools.chromsizes_to_file`
will create a file from a dictionary or string:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> # with no filename specified, a tempfile will be created
    >>> pybedtools.chromsizes_to_file(pybedtools.chromsizes('dm3'), 'dm3.genome')
    'dm3.genome'
    >>> print(open('dm3.genome').read())
    chr2L	23011544
    chr2R	21146708
    chr3L	24543557
    chr3R	27905053
    chr4	1351857
    chrX	22422827
    chr2LHet	368872
    chr2RHet	3288761
    chr3LHet	2555491
    chr3RHet	2517507
    chrM	19517
    chrU	10049037
    chrUextra	29004656
    chrXHet	204112
    chrYHet	347038
    <BLANKLINE>
