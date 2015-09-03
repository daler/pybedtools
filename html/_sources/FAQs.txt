.. include:: includeme.rst

FAQs
====

.. note::

    More detailed answers to these questions can often be found on the `Issues
    <https://github.com/daler/pybedtools/issues/>`_ page.

"Does pybedtools have a simple reader/writer for BED files?"
------------------------------------------------------------
While `pybedtools` designed to work with BEDTools, the reading/writing/parsing
function can be easily used for other things.

Simply iterating over a :class:`BedTool` object will parse each line into
a :class:`Interval` object.  You can then manipulate this or access the fields
as needed.

For example::

    x = pybedtools.example_bedtool('a.bed')
    for interval in x:
        # do something with interval

However, if you're planning on writing the results out to file, it may be more
useful to write a transformation function along with the :meth:`BedTool.each`
method.  This allows you to read, transform, and write all in one command::

    def my_func(f):
        """
        adds 10 bp to the stop
        """
        f.stop += 1
        return f

    pybedtools.example_bedtool('a.bed')\
        .each(my_func)\
        .saveas('out.bed')

Another useful idiom is creating a generator function.  For example, here we
change the name field to reflect the value of a counter.  We create a BedTool
from the iterator and then save it::

    def gen():
        counter = 0
        for i in pybedtools.example_bedtool('a.bed'):
            i.name = str(counter)
            counter += 1
            yield i

    pybedtools.BedTool(gen()).saveas('counted.bed')


See :ref:`saveresults` for more on saving the results.

"Can I create a BedTool object from an existing list?"
------------------------------------------------------

Sure, the :class:`BedTool` constructor will figure it out::

    items = [
        ('chr1', 100, 200),
        ('chr1', 500, 600),
    ]

    x = pybedtools.BedTool(items)


"I'm getting an empty BedTool"
------------------------------
Check to make sure you're not consuming a BedTool generator.  Note that
:meth:`BedTool.filter` and :meth:`BedTool.each` will return a generator BedTool
object. Keep in mind that checking the length of a generator BedTool will
completely consume it.

It's probably best to save intermediate versions to file using
:meth:`BedTool.saveas`.  If you don't provide a filename, it'll save to an
automatically cleaned up tempfile::

    my_bedtool\
     .filter(my_filter_func)\
     .saveas()\
     .intersect(y)\
     .filter(lambda x: len(x) > 1000)\
     .saveas('filtered-intersected-large.bed')


"I'm getting a MalformedBedLineError"
-------------------------------------
This error can be raised by BEDTools itself.  Typical reasons are that start
> end, or the fields are not tab-delimited.

You can try the :func:`pybedtools.remove_invalid` function to clean up your
file, or manually edit the offending lines.


"I get a segfault when iterating over a BedTool object"
-------------------------------------------------------

`Issue #88 <https://github.com/daler/pybedtools/issues/88>`_ which
addresses this issue -- in summary, Cython's handling of iterators works
unexpectedly.  It's best to call the `next()` method explicitly when doing
complex manipulations on an iterating :class:`BedTool`.


"Can I add extra information to FASTA headers when using BedTool.sequence()?"
-----------------------------------------------------------------------------

Since BEDTools adds the feature name to the FASTA header, you can manipulate
the feature name on the fly with a custom modifier function::

    def fields2name(f):
        "replace GFF featuretype field with the attributes field"
        f[2] = f[-1]
        return f

    import pybedtools
    g = pybedtools.BedTool("my.gff").each(fields2name).sequence(fi='my.fasta')

    print open(g.seqfn).readline()


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
