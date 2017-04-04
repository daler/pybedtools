Scripts
=======
:mod:`pybedtools` comes with several scripts that illustrate common use-cases.


In Python 2.7, you can use::

    python -m pybedtools

to get a list of scripts and their description.

Venn diagram scripts
--------------------
There are two scripts for making Venn diagrams, depending on how you'd like the
diagrams to look.  Both simply take 3 BED files as input.  ``venn_gchart.py``
uses the Google Chart API, while ``venn_mpl.py`` uses matplotlib if you have it
installed.

Upon installing :mod:`pybedtools`, these scripts should be available on your
path.  Calling them with the `-h` option will print the help, and using the
`--test` option will run a test, creating a new file `out.png` in the current
working directory.

.. figure:: images/gchart.png
    :width: 300px

    Above: using `--test` with `venn_gchart.py` results in this figure


.. figure:: images/mpl.png
    :width: 500px

    Above: Result of using `--test` with `venn_mpl.py`


Intron/exon classification
--------------------------
The script `intron_exon_reads.py` accepts a GFF file (with introns and exons
annotated) and a BAM file.  When complete, it prints out the number of exonic,
intronic, and both intronic and exonic (i.e., from overlapping genes or
isoforms).  This script is also a good example of how to do use Python's
:mod:`multiprocessing` for parallel computation.

Annotate.py
-----------
The `annotate.py` script extends `closestBed` by classifying features (intron,
exon) that are a distance of 0 away from the query features.
