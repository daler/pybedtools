import os
import sys
import itertools
import six
import time
import pysam
import pybedtools


def tag_bedpe(bedpe, queries, verbose=False):
    """
    Tag each end of a BEDPE with a set of (possibly many) query BED files.

    For example, given a BEDPE of interacting fragments from a Hi-C experiment,
    identify the contacts between promoters and ChIP-seq peaks. In this case,
    promoters and ChIP-seq peaks of interest would be provided as BED files.

    The strategy is to split the BEDPE into two separate files.  Each file is
    intersected independently with the set of queries.  The results are then
    iterated through in parallel to tie the ends back together. It is this
    iterator that is returned (see example below).

    Parameters
    ----------

    bedpe : str
        BEDPE-format file. Must be name-sorted.

    queries : dict
        Dictionary of BED/GFF/GTF/VCF files to use. After splitting the BEDPE,
        these query files (values in the dictionary) will be passed as the `-b`
        arg to `bedtools intersect`. The keys are passed as the `names`
        argument for `bedtools intersect`

        *Features in each file must have unique names*. Use
        :func:`pybedtools.featurefuncs.UniqueID` to help fix this.

        *Each file must be BED3 to BED6*.

    Returns
    -------
    Tuple of (iterator, n, extra).

    `iterator` is described below. `n` is the total number of lines in the
    BEDPE file, which is useful for calculating percentage complete for
    downstream work. `extra` is the number of extra fields found in the BEDPE
    (also useful for downstream processing).

    `iterator` yields tuples of (label, end1_hits, end2_hits) where `label` is
    the name field of one line of the original BEDPE file. `end1_hits` and
    `end2_hits` are each iterators of BED-like lines representing all
    identified intersections across all query BED files for end1 and end2 for
    this pair.

    Recall that BEDPE format defines a single name and a single score for each
    pair. For each item in `end1_hits`, the fields are::

        chrom1
        start1
        end1
        name
        score
        strand1
        [extra fields]
        query_label
        fields_from_query_intersecting_end1

    where `[extra fields]` are any additional fields from the original BEDPE,
    `query_label` is one of the keys in the `beds` input dictionary, and the
    remaining fields in the line are the intersecting line from the
    corresponding BED file in the `beds` input dictionary.

    Similarly, each item in `end2_hits` consists of:

        chrom2
        start2
        end2
        name
        score
        strand2
        [extra fields]
        query_label
        fields_from_query_intersecting_end2

    At least one line is reported for every line in the BEDPE file. If there
    was no intersection, the standard BEDTools null fields will be shown. In
    `end1_hits` and `end2_hits`, a line will be reported for each hit in each
    query.

    Example
    -------

    Consider the following BEDPE (where "x1" is an aribtrary extra field).

    >>> bedpe = pybedtools.example_bedtool('test_bedpe.bed')
    >>> print(bedpe) # doctest: +NORMALIZE_WHITESPACE
    chr1  1  10  chr1  50   90   pair1  5  +  -  x1
    chr1  2  15  chr1  200  210  pair2  1  +  +  y1
    <BLANKLINE>


    And the following transcription start sites (TSSes) in BED4 format:

    >>> tsses = pybedtools.example_bedtool('test_tsses.bed')
    >>> print(tsses) # doctest: +NORMALIZE_WHITESPACE
    chr1  5   6   gene1
    chr1  60  61  gene2
    chr1  88  89  gene3
    <BLANKLINE>

    And the following called peaks as BED6:

    >>> peaks = pybedtools.example_bedtool('test_peaks.bed')
    >>> print(peaks) # doctest: +NORMALIZE_WHITESPACE
    chr1  3  4  peak1  50  .
    <BLANKLINE>

    Then we can get the following iterator, n, and extra. Note that the
    OrderedDict is only for testing to ensure output is always consistend; in
    practice a regular dictionary is fine:

    >>> from pybedtools.contrib.long_range_interaction import tag_bedpe
    >>> from collections import OrderedDict
    >>> queries = OrderedDict()
    >>> queries['tss'] = tsses
    >>> queries['pk'] = peaks
    >>> iterator, n, extra = tag_bedpe(bedpe, queries)
    >>> print(n)
    2
    >>> print(extra)
    1

    The following illustrates that each item in the iterator represents one
    pair, and each item in each group represents an intersection with one end.
    Note that the sorting is necessary only for the doctests to be output in
    consistent format; this not typically needed:

    >>> for (label, end1_hits, end2_hits) in iterator:
    ...    end1_hits = sorted(end1_hits, key=lambda x: str(x))
    ...    end2_hits = sorted(end2_hits, key=lambda x: str(x))
    ...    print('PAIR = {}'.format(label))
    ...    print('end1_hits:')
    ...    for i in end1_hits:
    ...        print(i, end='')
    ...    print('end2_hits:')
    ...    for i in end2_hits:
    ...        print(i, end='')  # doctest: +NORMALIZE_WHITESPACE
    PAIR = pair1
    end1_hits:
    chr1       1       10      pair1   5       +       x1      pk      chr1    3       4       peak1   50      .       1
    chr1       1       10      pair1   5       +       x1      tss     chr1    5       6       gene1   1
    end2_hits:
    chr1       50      90      pair1   5       -       x1      tss     chr1    60      61      gene2   1
    chr1       50      90      pair1   5       -       x1      tss     chr1    88      89      gene3   1
    PAIR = pair2
    end1_hits:
    chr1       2       15      pair2   1       +       y1      pk      chr1    3       4       peak1   50      .       1
    chr1       2       15      pair2   1       +       y1      tss     chr1    5       6       gene1   1
    end2_hits:
    chr1       200     210     pair2   1       +       y1      .       .       -1      -1      .       0

    See the `cis_trans_interactions()` function for one way of summarizing
    these data.
    """
    b = pybedtools.BedTool(bedpe)

    # Figure out if the supplied bedpe had any extra fields. If so, the fields
    # are repeated in each of the split output files.
    observed = b.field_count()
    extra = observed - 10
    extra_inds = [10 + i for i in range(extra)]

    end1_fn = pybedtools.BedTool._tmp()
    end2_fn = pybedtools.BedTool._tmp()

    # Performance notes:
    # We don't need the overhead of converting every line into
    # a pybedtools.Interval object just so we can grab the fields. Doing so
    # takes 3.5x more time than simply splitting each line on a tab.
    if verbose:
        print("splitting BEDPE into separate files.")
        print("end1 is going to %s" % end1_fn)
        print("end2 is going to %s" % end2_fn)

    n = 0
    with open(end1_fn, "w") as end1_out, open(end2_fn, "w") as end2_out:
        for line in open(b.fn):
            n += 1
            f = line.strip().split("\t")
            end1_out.write(
                "\t".join((f[i] for i in [0, 1, 2, 6, 7, 8] + extra_inds)) + "\n"
            )
            end2_out.write(
                "\t".join((f[i] for i in [3, 4, 5, 6, 7, 9] + extra_inds)) + "\n"
            )

    # Performance notes:
    #
    # For small BEDPE and large set of query files, it would be faster to sort
    # these independently, intersect with sorted=True, and then re-sort by name
    # for the grouping. For large BEDPE, I don't think the sorted=True
    # performance gain outweighs the hit from sorting twice.
    #
    # On the other hand, if BEDPE was coord-sorted in the first place, only
    # end2 would need to be sorted and re-sorted. On the other (third!?) hand,
    # BEDPE creation from BAM implies name-sorting, so it's probably not
    # reasonable to assume coord-sorted.
    #
    # In the end: don't do any sorting.

    end1_bt = pybedtools.BedTool(end1_fn)
    end2_bt = pybedtools.BedTool(end2_fn)
    names, fns = [], []
    for name, fn in queries.items():
        names.append(name)
        if isinstance(fn, pybedtools.BedTool):
            fns.append(fn.fn)
        else:
            fns.append(fn)

    if verbose:
        print("intersecting end 1")
    end1_hits = end1_bt.intersect(list(fns), names=names, wao=True)
    if verbose:
        print("intersecting end 2")
    end2_hits = end2_bt.intersect(list(fns), names=names, wao=True)
    if verbose:
        print("intersection with end1 is in %s" % (end1_hits.fn))
        print("intersection with end2 is in %s" % (end2_hits.fn))

    grouped_end1 = itertools.groupby(end1_hits, lambda f: f[3])
    grouped_end2 = itertools.groupby(end2_hits, lambda f: f[3])

    def gen():
        for (label1, group1), (label2, group2) in six.moves.zip(
            grouped_end1, grouped_end2
        ):
            assert label1 == label2
            yield label1, group1, group2

    return gen(), n, extra


def cis_trans_interactions(iterator, n, extra, verbose=True):
    """
    Converts the output from `tag_bedpe` into a pandas DataFrame containing
    information about regions that contact each other in cis (same fragment) or
    trans (different fragments).

    For example, given a BEDPE file representing 3D interactions in the genome,
    we want to identify which transcription start sites are connected to distal
    regions containing a peak.

    >>> bedpe = pybedtools.example_bedtool('test_bedpe.bed')
    >>> print(bedpe) # doctest: +NORMALIZE_WHITESPACE
    chr1  1  10  chr1  50   90   pair1  5  +  -  x1
    chr1  2  15  chr1  200  210  pair2  1  +  +  y1
    <BLANKLINE>

    >>> tsses = pybedtools.example_bedtool('test_tsses.bed')
    >>> print(tsses) # doctest: +NORMALIZE_WHITESPACE
    chr1  5   6   gene1
    chr1  60  61  gene2
    chr1  88  89  gene3
    <BLANKLINE>

    >>> peaks = pybedtools.example_bedtool('test_peaks.bed')
    >>> print(peaks) # doctest: +NORMALIZE_WHITESPACE
    chr1  3  4  peak1  50  .
    <BLANKLINE>

    Here's what the tracks look like. Note that pair1 is evidence of
    a gene1-gene2 interaction and a gene1-gene3 interaction::

        TRACKS:

                      1         2  /         5         6    / 8         9     / 20
            0123456789012345678901 / 2345678901234567890123 / 012345678901234 / 0123456789
       pair1 |||||||||------------ / --------|||||||||||||| / ||||||||||||||| /
       pair2  |||||||||||||------- / ---------------------- / --------------- / ||||||||||
       tsses     1                 /                   2    /         3
       peaks   1


    >>> from collections import OrderedDict
    >>> queries = OrderedDict()
    >>> queries['tss'] = tsses
    >>> queries['pk'] = peaks
    >>> iterator, n, extra = tag_bedpe(bedpe, queries)
    >>> for (label, group1, group2) in iterator:
    ...    group1 = sorted(group1, key=lambda x: str(x))
    ...    group2 = sorted(group2, key=lambda x: str(x))
    ...    for i in group1:
    ...        print(i, end='')  # doctest: +NORMALIZE_WHITESPACE
    ...    for i in group2:
    ...        print(i, end='')  # doctest: +NORMALIZE_WHITESPACE
    chr1       1       10      pair1   5       +       x1      pk      chr1    3       4       peak1   50      .       1
    chr1       1       10      pair1   5       +       x1      tss     chr1    5       6       gene1   1
    chr1       50      90      pair1   5       -       x1      tss     chr1    60      61      gene2   1
    chr1       50      90      pair1   5       -       x1      tss     chr1    88      89      gene3   1
    chr1       2       15      pair2   1       +       y1      pk      chr1    3       4       peak1   50      .       1
    chr1       2       15      pair2   1       +       y1      tss     chr1    5       6       gene1   1
    chr1       200     210     pair2   1       +       y1      .       .       -1      -1      .       0

    Now we run the same thing, but now aggregate it. Note that each piece of
    interaction evidence has its own line. The first line shows that pair1 has
    gene1 and peak1 in the same fragment, and that they are connected to gene2.
    The second line shows again that gene1 and peak1 are in the same fragmet
    and that they are also connected to gene3:

    >>> import pandas; pandas.set_option('display.max_columns', 10)
    >>> iterator, n, extra = tag_bedpe(bedpe, {'tss': tsses, 'pk': peaks})
    >>> df =  cis_trans_interactions(iterator, n, extra)
    >>> print(df.sort_values(list(df.columns)).reset_index(drop=True))
      target_label target_name cis_label cis_name distal_label distal_name  label
    0           pk       peak1       tss    gene1            .           .  pair2
    1           pk       peak1       tss    gene1          tss       gene2  pair1
    2           pk       peak1       tss    gene1          tss       gene3  pair1
    3          tss       gene1        pk    peak1            .           .  pair2
    4          tss       gene1        pk    peak1          tss       gene2  pair1
    5          tss       gene1        pk    peak1          tss       gene3  pair1
    6          tss       gene2       tss    gene3           pk       peak1  pair1
    7          tss       gene2       tss    gene3          tss       gene1  pair1
    8          tss       gene3       tss    gene2           pk       peak1  pair1
    9          tss       gene3       tss    gene2          tss       gene1  pair1

    If we only care about genes:

    >>> print((df[df.target_label == 'tss']).sort_values(list(df.columns)).reset_index(drop=True))
      target_label target_name cis_label cis_name distal_label distal_name  label
    0          tss       gene1        pk    peak1            .           .  pair2
    1          tss       gene1        pk    peak1          tss       gene2  pair1
    2          tss       gene1        pk    peak1          tss       gene3  pair1
    3          tss       gene2       tss    gene3           pk       peak1  pair1
    4          tss       gene2       tss    gene3          tss       gene1  pair1
    5          tss       gene3       tss    gene2           pk       peak1  pair1
    6          tss       gene3       tss    gene2          tss       gene1  pair1


    Note that in pair2, there is no evidence of interaction between gene1 and
    gene2.

    What interacts distally with gene2's TSS?

    >>> assert set(df.loc[df.target_name == 'gene2', 'distal_name']).difference('.') == set([u'gene1', u'peak1'])

    """
    try:
        import pandas
    except ImportError:
        raise ImportError("pandas must be installed to use this function")
    c = 0
    lines = []
    for label, end1_hits, end2_hits in iterator:
        c += 1
        if c % 1000 == 0:
            print("%d (%.1f%%)\r" % (c, c / float(n) * 100), end="")
            sys.stdout.flush()

        # end1_hits has the full lines of all intersections with end1
        end1_hits = list(end1_hits)
        end2_hits = list(end2_hits)

        def get_name_hits(f):
            """
            Returns the key (from which file the interval came) and the name
            (of the individual feature).
            """
            # this is the "name" reported if there was no hit.
            if f[6 + extra] == ".":
                return (".", ".")
            interval = pybedtools.create_interval_from_list(f[7 + extra :])
            return [f[6 + extra], interval.name]

        names1 = set(six.moves.map(tuple, six.moves.map(get_name_hits, end1_hits)))
        names2 = set(six.moves.map(tuple, six.moves.map(get_name_hits, end2_hits)))

        for cis, others in [(names1, names2), (names2, names1)]:
            for target in cis:
                if target == (".", "."):
                    continue
                non_targets = set(cis).difference([target])
                if len(non_targets) == 0:
                    non_targets = [(".", ".")]
                for non_target in non_targets:
                    for other in others:
                        line = []
                        line.extend(target)
                        line.extend(non_target)
                        line.extend(other)
                        line.append(label)
                        lines.append(line)

    df = pandas.DataFrame(
        lines,
        columns=[
            "target_label",
            "target_name",
            "cis_label",
            "cis_name",
            "distal_label",
            "distal_name",
            "label",
        ],
    )
    df = df.drop_duplicates()
    return df
