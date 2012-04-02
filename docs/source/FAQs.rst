.. include:: includeme.rst

FAQs
====

"Too many files open" error
---------------------------

Sometimes you may get the error::

    * Too many files open -- please submit a bug report so that this can be fixed

This error occurs because you have hit your operating system's limit on the
number of open files.  This usually happens when creating many :class:`BedTool`
objects, often within a for-loop.

In general, **try to create as few** :class:`BedTool` **objects as you can**.  Every time you
create a :class:`BedTool` object, you create a new open file.  There is usually
a BEDTools program that already does what you want, and will do it faster.


For example, say we want to:

* start with all annotations
* only consider exons
* write a file containing just exons
* count reads in multiple BAM files for each exon


Here is a first draft.  Note that the for-loop creates a :class:`BedTool`
object each iteration, and the `result` is yet another :class:`BedTool`.  This
will version will raise the "Too many files open" error.

.. code-block:: python

    # This version will be slow and, with many exons, will raise the "Too many
    # files open" error

    import pybedtools
    all_features = pybedtools.BedTool('annotations.gff')
    fout = open('exons.gff', 'w')
    for feature in all_features:
        if feature[2] != 'exon':
            continue

        fout.write(str(feature))

        bt = pybedtools.BedTool([feature])
        result = bt.multi_bam_coverage(bams=['reads1.bam', 'reads2.bam'])

        # ...do something with result

    fout.close()

In contrast, it would be better to construct an "exon-only" :class:`BedTool` at
the beginning.  The :meth:`BedTool.filter` method is a good way to do this.
Then, there is only one call to :meth:`BedTool.multi_bam_coverage`.

In this version there are only 3 :class:`BedTool` objects:  the
one that opens `annotations.gff`, the one that uses `exons.gff` after it is
saved, and `result`.  (Note that the one created from the filter operation is
a "streaming" BedTool, so there is no open file that will contribute to the
total).

.. code-block:: python

    # This is the recommended way.

    import pybedtools

    exons = pybedtools.BedTool('annotations.gff')\
        .filter(lambda x: x[2] == 'exon')\
        .saveas('exons.gff')

    result = exons.multi_bam_coverage(bams=['reads1.bam', 'reads2.bam'])

    # ...do something with result
