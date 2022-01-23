import pybedtools
import gzip
import os
import subprocess
import sys
from textwrap import dedent
import six
import pytest
import psutil


testdir = os.path.dirname(__file__)
tempdir = os.path.join(os.path.abspath(testdir), "tmp")
unwriteable = "unwriteable"


def setup_module():
    if not os.path.exists(tempdir):
        os.system("mkdir -p %s" % tempdir)
    pybedtools.set_tempdir(tempdir)


def teardown_module():
    if os.path.exists(tempdir):
        os.system("rm -r %s" % tempdir)
    pybedtools.cleanup()


def fix(x):
    """
    Replaces spaces with tabs, removes spurious newlines, and lstrip()s each
    line. Makes it really easy to create BED files on the fly for testing and
    checking.
    """
    s = ""
    for i in x.splitlines():
        i = i.lstrip()
        if i.endswith("\t"):
            add_tab = "\t"
        else:
            add_tab = ""
        if len(i) == 0:
            continue
        i = i.split()
        i = "\t".join(i) + add_tab + "\n"
        s += i
    return s


def test_issue_81():
    genome = {"chr1": (0, 5000)}
    result = pybedtools.BedTool().window_maker(genome=genome, w=1000, s=500)
    assert result == fix(
        """
        chr1	0	1000
        chr1	500	1500
        chr1	1000	2000
        chr1	1500	2500
        chr1	2000	3000
        chr1	2500	3500
        chr1	3000	4000
        chr1	3500	4500
        chr1	4000	5000
        chr1	4500	5000
        """
    ), result


def test_issue_118():
    p = psutil.Process(os.getpid())
    start_fds = p.num_fds()
    a = pybedtools.example_bedtool("a.bed")
    b = pybedtools.example_bedtool("b.bed")
    for i in range(100):
        c = a.intersect(b)
        c.field_count()
    stop_fds = p.num_fds()
    assert start_fds == stop_fds


def test_issue_131():
    """
    Regression test; in previous versions this would cause a segfault.
    """
    from itertools import groupby

    x = pybedtools.BedTool(
        [
            ("chr1", 12, 13, "N", 1000, "+"),
            ("chr1", 12, 13, "N", 1000, "-"),
            ("chr1", 12, 13, "N", 1000, "-"),
            ("chr1", 115, 116, "N", 1000, "+"),
        ]
    )

    for key, group_ in groupby(x, key=lambda r: (r.chrom, r.start, r.end)):
        print(key, map(lambda r: r.strand, group_))


def test_issue_138():
    x = pybedtools.BedTool(
        """
        chr1	1	100
        """,
        from_string=True,
    )
    y = pybedtools.BedTool(
        """
        chr2	500	600	feature1
        """,
        from_string=True,
    )
    z = pybedtools.BedTool(
        """
        chr3	99	999	feature2	0	+
        """,
        from_string=True,
    )

    # When force_truncate is False (default), output field length should be the
    # min across all input field lengths.
    assert y.cat(y, z, postmerge=False) == fix(
        """
        chr2	500	600	feature1
        chr2	500	600	feature1
        chr3	99	999	feature2
        """
    )

    assert y.cat(x, z, postmerge=False) == fix(
        """
        chr2	500	600
        chr1	1	100
        chr3	99	999
        """
    )

    assert z.cat(z, z, z, z, postmerge=False) == fix(
        """
        chr3	99	999	feature2	0	+
        chr3	99	999	feature2	0	+
        chr3	99	999	feature2	0	+
        chr3	99	999	feature2	0	+
        chr3	99	999	feature2	0	+
        """
    )


def test_issue_141():
    a = pybedtools.example_bedtool("a.bed")
    b = pybedtools.example_bedtool("b.bed")

    # make an empty file
    empty = pybedtools.BedTool("", from_string=True)

    # invalid file format
    malformed = pybedtools.BedTool("a	a	a", from_string=True)

    # positive control; works
    a + b

    # "adding" an empty file always gets zero features
    assert len(a + empty) == 0
    assert len(empty + a) == 0
    assert len(empty + empty) == 0

    # "adding" a malformed file raises MalformedBedLineError
    # (an uncaught exception raised when trying to intersect)
    with pytest.raises(pybedtools.MalformedBedLineError):
        a + malformed

    x = pybedtools.example_bedtool("x.bam")
    x + a


def test_issue_141b():
    a = pybedtools.example_bedtool("hg38-problem.bed")
    b = pybedtools.example_bedtool("hg38-base.bed")

    # prior to fixing #147, BEDToolsError was raised here due to unhandled
    # stderr.  Now the stderr is detected as OK because it's just a warning, so
    # these lines are commented out now.
    # assert_raises(pybedtools.helpers.BEDToolsError, a.intersect, b)
    # assert_raises(pybedtools.helpers.BEDToolsError, a.__add__, b)

    # use nonamecheck
    res = a.intersect(b, nonamecheck=True)
    assert res == fix(
        """
        chr1 2 50
        """
    )


def test_issue_143():
    def func(x):
        x.start += 10
        return x

    a = pybedtools.example_bedtool("a.bed")
    b = a.merge(s=True, stream=True).each(func).saveas()
    c = a.merge(s=True).each(func).saveas()
    assert b == c

    b = a.merge(s=True, stream=True)
    for i in b:
        assert isinstance(i, pybedtools.Interval)

    b = a.merge(s=True, stream=True)
    for i in iter(iter(iter(b))):
        assert isinstance(i, pybedtools.Interval)

    for i in a.merge(s=True, stream=True).each(lambda x: x):
        assert isinstance(i, pybedtools.Interval)


def test_issue_145():
    x = pybedtools.BedTool(
        """
    chr1    1   100 feature1    0   +
    chr1    1   100 feature1    0   +
    """,
        from_string=True,
    ).saveas("foo.bed")

    g = pybedtools.chromsizes_to_file({"chr1": (0, 200)}, "genome.txt")
    y = x.genome_coverage(g=g, **{"5": True})

    # trying to print causes pybedtools to interpret as a BED file, but it's
    # a histogram so line 2 raises error
    with pytest.raises(pybedtools.MalformedBedLineError):
        print(y)

    # solution is to iterate over lines of file; make sure this works
    for line in open(y.fn):
        print(line)

    # if streaming, iterate over y.fn directly:
    y = x.genome_coverage(g=g, **{"5": True})
    for line in y.fn:
        print(line)


def test_issue_147():
    # previously this would raise BEDToolsError because of unexpected stderr.
    with open(pybedtools.BedTool._tmp(), "w") as tmp:
        orig_stderr = sys.stderr
        sys.stderr = tmp
        v = pybedtools.example_bedtool("vcf-stderr-test.vcf")
        b = pybedtools.example_bedtool("vcf-stderr-test.bed")
        v.intersect(b)
    sys.stderr = orig_stderr


def test_issue_154():
    regions = [("chr2", int(1), int(2), "tag")]
    pybedtools.BedTool(regions)

    # ensure longs are OK as start/stop. In Python3 everything is long, so only
    # try the following on PY2
    if six.PY2:
        regions = [("chr2", long(1), long(2), "tag")]
        pybedtools.BedTool(regions)


def test_issue_151():
    # this used to be incorrectly inferred to be SAM because of the name field
    # being an integer and >=11 fields. The fix was to check the strand -- if
    # it's not in ['+', '-', '.'] then consider it a SAM.
    f = pybedtools.create_interval_from_list(
        [
            "chr1",
            "1197700",
            "1197758",
            "0",
            "0.318355266754715",
            "-",
            "foo",
            "bar",
            "bam",
            "baz",
            "bizzle",
            "buz",
            "bis",
        ]
    )
    assert f.file_type == "bed"


def test_issue_156():
    # NOTE: this isn't appropriate for including in the test_iter cases, since
    # that tests filenames, gzipped files, and iterators. There's no support
    # for "list of iterators" as the `b` argument. Plus, here we're not
    # concerned with the ability to handle those different input types -- just
    # that lists of filenames works.
    a = pybedtools.example_bedtool("a.bed")
    b = [pybedtools.example_filename("b.bed"), pybedtools.example_filename("c.gff")]
    res = str(a.intersect(b))
    assert res == fix(
        """
        chr1    59      100     feature1        0       +
        chr1    155     200     feature2        0       +
        chr1    173     200     feature2        0       +
        chr1    173     200     feature2        0       +
        chr1    100     200     feature2        0       +
        chr1    155     200     feature3        0       -
        chr1    464     500     feature3        0       -
        chr1    485     500     feature3        0       -
        chr1    173     326     feature3        0       -
        chr1    438     500     feature3        0       -
        chr1    495     500     feature3        0       -
        chr1    485     500     feature3        0       -
        chr1    173     326     feature3        0       -
        chr1    438     500     feature3        0       -
        chr1    150     269     feature3        0       -
        chr1    900     901     feature4        0       +
        chr1    900     913     feature4        0       +
        chr1    900     913     feature4        0       +
        chr1    900     950     feature4        0       +
        """
    ), res

    res = str(a.intersect(b, wb=True, names=["B", "C"]))
    assert res == fix(
        """
        chr1	59	100	feature1	0	+	C	chr1	ucb	gene	60	269	.	-	.	ID=thaliana_1_6160_6269;match=fgenesh1_pg.C_scaffold_1000119;rname=thaliana_1_6160_6269
        chr1	155	200	feature2	0	+	B	chr1	155	200	feature5	0	-
        chr1	173	200	feature2	0	+	C	chr1	ucb	CDS	174	326	.	+	.	Parent=AT1G01010.mRNA;rname=AT1G01010
        chr1	173	200	feature2	0	+	C	chr1	ucb	mRNA	174	326	.	+	.	ID=AT1G01010.mRNA;Parent=AT1G01010;rname=AT1G01010
        chr1	100	200	feature2	0	+	C	chr1	ucb	gene	60	269	.	-	.	ID=thaliana_1_6160_6269;match=fgenesh1_pg.C_scaffold_1000119;rname=thaliana_1_6160_6269
        chr1	155	200	feature3	0	-	B	chr1	155	200	feature5	0	-
        chr1	464	500	feature3	0	-	C	chr1	ucb	gene	465	805	.	+	.	ID=thaliana_1_465_805;match=scaffold_801404.1;rname=thaliana_1_465_805
        chr1	485	500	feature3	0	-	C	chr1	ucb	CDS	486	605	.	+	.	Parent=AT1G01010.mRNA;rname=AT1G01010
        chr1	173	326	feature3	0	-	C	chr1	ucb	CDS	174	326	.	+	.	Parent=AT1G01010.mRNA;rname=AT1G01010
        chr1	438	500	feature3	0	-	C	chr1	ucb	CDS	439	630	.	+	.	Parent=AT1G01010.mRNA;rname=AT1G01010
        chr1	495	500	feature3	0	-	C	chr1	ucb	mRNA	496	576	.	+	.	ID=AT1G01010.mRNA;Parent=AT1G01010;rname=AT1G01010
        chr1	485	500	feature3	0	-	C	chr1	ucb	mRNA	486	605	.	+	.	ID=AT1G01010.mRNA;Parent=AT1G01010;rname=AT1G01010
        chr1	173	326	feature3	0	-	C	chr1	ucb	mRNA	174	326	.	+	.	ID=AT1G01010.mRNA;Parent=AT1G01010;rname=AT1G01010
        chr1	438	500	feature3	0	-	C	chr1	ucb	mRNA	439	899	.	+	.	ID=AT1G01010.mRNA;Parent=AT1G01010;rname=AT1G01010
        chr1	150	269	feature3	0	-	C	chr1	ucb	gene	60	269	.	-	.	ID=thaliana_1_6160_6269;match=fgenesh1_pg.C_scaffold_1000119;rname=thaliana_1_6160_6269
        chr1	900	901	feature4	0	+	B	chr1	800	901	feature6	0	+
        chr1	900	913	feature4	0	+	C	chr1	ucb	mRNA	631	913	.	+	.	ID=AT1G01010.mRNA;Parent=AT1G01010;rname=AT1G01010
        chr1	900	913	feature4	0	+	C	chr1	ucb	CDS	760	913	.	+	.	Parent=AT1G01010.mRNA;rname=AT1G01010
        chr1	900	950	feature4	0	+	C	chr1	ucb	CDS	706	1095	.	+	.	Parent=AT1G01010.mRNA;rname=AT1G01010
        """
    ), res


def test_issue_157():
    # the problem here was that converting to file from dataframe didn't pass
    # through enough options to pandas.
    try:
        import pandas
    except ImportError:
        pytest.xfail("pandas not installed; skipping test")
    vcf = pybedtools.example_bedtool("1000genomes-example.vcf")
    bed = pybedtools.BedTool("20\t14300\t17000", from_string=True)
    non_dataframe = str(vcf.intersect(bed))
    df = vcf.to_dataframe(
        comment="#",
        names=[
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "NA00001",
            "NA00002",
            "NA00003",
        ],
    )

    header = "".join([line for line in open(vcf.fn) if line.startswith("#")])
    outfile = pybedtools.BedTool._tmp()
    with open(outfile, "w") as fout:
        fout.write(header)
        vcf_from_df = pybedtools.BedTool.from_dataframe(df, outfile=fout)
    from_dataframe = str(vcf_from_df.intersect(bed))
    assert non_dataframe == from_dataframe


def test_PR_158():
    # See #121 for original, #122 for follow-up, and #158 for fix.
    #
    # This used to crash with "OverflowError: can't convert negative value to CHRPOS"
    b = pybedtools.example_bedtool("issue_121.bam")
    print(b)


def test_issue_162():
    a = pybedtools.BedTool("", from_string=True)
    b = pybedtools.example_bedtool("b.bed")
    c = pybedtools.BedTool()
    with pytest.raises(ValueError):
        b.cat(c)
    assert str(b.cat(a)) == fix(
        """
        chr1	155	200
        chr1	800	901
        """
    )


def test_issue_164():
    a = pybedtools.example_bedtool("164.gtf")
    y = a.filter(
        lambda gene: gene.name in ["ENSMUSG00000000003", "ENSMUSG00000000037"]
    ).saveas()
    # don't use the fix() convenience function because we have both tabs (field
    # sep) and spaces (attributes sep)
    expected = dedent(
        """\
    chrX	gffutils_derived	gene	77837901	77853623	.	-	.	gene_id "ENSMUSG00000000003";
    chrX	gffutils_derived	gene	161117193	161258213	.	+	.	gene_id "ENSMUSG00000000037";
    """
    )
    assert str(y) == expected


def test_issue_168():
    # Regression test:
    # this would previously segfault in at least pysam 0.8.4
    #
    x = pybedtools.example_bedtool("1000genomes-example.vcf")
    fn = x.bgzip(is_sorted=True, force=True)
    y = pybedtools.BedTool(fn)


def test_issue_169():
    x = pybedtools.example_bedtool("1000genomes-example.vcf")
    fn = x.bgzip(is_sorted=False, force=True)
    line = gzip.open(fn, "rt").readline()
    assert str(line).startswith("#"), line


def test_issue_196():
    bed = pybedtools.BedTool(
        """
        8 129185980 129186130 A 0.1
        8 129185980 129186130 B 0.2
        """,
        from_string=True,
    )
    bed = bed.tabix()
    snp = pybedtools.BedTool("8\t129186110\t129186111\trs72722756", from_string=True)
    intersection = bed.tabix_intervals(
        "{}:{}-{}".format("8", 129186110, 129186111)
    ).intersect(snp, wa=True, wb=True)

    # prior to fixing this issue, intervals would be concatenated. This was
    # because pysam.ctabix.tabixIterator does not include newlines when
    # yielding. The incorrect output was this:
    """
    8   129185980   129186130   A   0.18   129185980   129186130   B   0.2   8   129186110   129186111   rs72722756
    """

    # but should be this:
    assert intersection == fix(
        """
        8       129185980       129186130       A       0.1     8       129186110       129186111       rs72722756
        8       129185980       129186130       B       0.2     8       129186110       129186111       rs72722756
        """
    )


def test_issue_178():
    try:
        fn = pybedtools.example_filename("gdc.othersort.bam")
        pybedtools.contrib.bigwig.bam_to_bigwig(fn, genome="dm3", output="tmp.bw")
        x = pybedtools.contrib.bigwig.bigwig_to_bedgraph("tmp.bw")
        assert x == fix(
            """
            chr2L   70      75      1
            chr2L   140     145     1
            chr2L   150     155     1
            chr2L   160     165     1
            chr2L   210     215     1
            chrX    10      15      1
            chrX    70      75      1
            chrX    140     145     1
            """
        )
        os.unlink("tmp.bw")

    # If bedGraphToBigWig is not on the path, see
    # https://github.com/daler/pybedtools/issues/227
    except FileNotFoundError:
        pass


def test_issue_180():
    a = pybedtools.example_bedtool("a.bed")
    a = a.tabix(force=True)
    assert a.tabix_contigs() == ["chr1"]


def test_issue_181():
    a = pybedtools.example_bedtool("a.bed")
    a = a.tabix(force=True)
    a.tabix_intervals("none:1-5")
    with pytest.raises(ValueError):
        a.tabix_intervals("none:1-5", check_coordinates=True)


def test_issue_203():
    x = pybedtools.example_bedtool("x.bed")
    x.truncate_to_chrom(genome="hg19")


def test_issue_217():

    # the doctest at
    # https://daler.github.io/pybedtools/intervals.html#common-interval-attributes
    # passes, so let's start with that
    x = pybedtools.example_bedtool("a.bed")[0]
    print(x)
    assert x.name == "feature1"

    # construct another interval using the same fields
    y = pybedtools.Interval(
        x.fields[0],
        int(x.fields[1]),
        int(x.fields[2]),
        x.fields[3],
        x.fields[4],
        x.fields[5],
    )

    # and using create_interval_from_list, which many internal functions use:
    z = pybedtools.create_interval_from_list(x.fields)

    # They are identical in all meaningful ways . . . .
    assert x.fields == y.fields == z.fields
    assert str(x) == str(y) == str(z) == "chr1\t1\t100\tfeature1\t0\t+\n"
    assert type(x) == type(y) == type(z) == pybedtools.Interval
    assert x.chrom == y.chrom == z.chrom == "chr1"
    assert x.start == y.start == z.start == 1
    assert x.stop == y.stop == z.stop == 100
    assert x.strand == y.strand == z.strand == "+"
    assert x.score == y.score == z.score == "0"

    assert x.file_type == "bed"

    # Previously this returned None
    assert y.file_type == "bed"

    assert z.file_type == "bed"

    assert x.name == "feature1"

    # Previously the directly-created Interval object returned None for a name.
    assert y.name == "feature1"

    assert z.name == "feature1"


def test_issue_218():
    from pybedtools.helpers import set_bedtools_path, get_bedtools_path
    from pybedtools import BedTool

    orig_path = get_bedtools_path()

    # As pointed out in #222, example_bedtool behaves differently from BedTool.
    # example_bedtool is defined in pybedtools.bedtool but pybedtools.BedTool
    # is imported in pybedtools.__init__. So check various constructors here.
    for constructor in (
        lambda x: pybedtools.example_bedtool(x),
        lambda x: pybedtools.BedTool(pybedtools.example_filename(x)),
        lambda x: pybedtools.bedtool.BedTool(pybedtools.example_filename(x)),
        # NOTE: we likely need recursive reloading (like IPython.deepreload)
        # for this to work:
        #
        # lambda x: BedTool(pybedtools.example_filename(x)),
    ):

        x = constructor("x.bed")
        x.sort()
        assert "Original BEDTools help" in pybedtools.bedtool.BedTool.sort.__doc__
        assert "Original BEDTools help" in x.sort.__doc__

        set_bedtools_path("nonexistent")

        # Calling BEDTools with non-existent path, but the docstring should not
        # have been changed.
        with pytest.raises(OSError):
            x.sort()
        assert "Original BEDTools help" in x.sort.__doc__

        # The class's docstring should have been reset though.
        assert pybedtools.bedtool.BedTool.sort.__doc__ is None

        # Creating a new BedTool object now that bedtools is not on the path
        # should detect that, adding a method that raises
        # NotImplementedError...
        y = constructor("x.bed")
        with pytest.raises(NotImplementedError):
            y.sort()

        # ...and correspondingly no docstring
        assert y.sort.__doc__ is None
        assert pybedtools.bedtool.BedTool.sort.__doc__ is None

        # Reset the path, and ensure the resetting works
        set_bedtools_path()
        z = constructor("x.bed")
        z.sort()


def test_issue_231():
    def filt(f):
        raise ValueError("failed")

    a = pybedtools.example_bedtool("a.bed")

    # Previously, ValueErrors in filter/each functions were silently ignored
    with pytest.raises(ValueError):
        assert list(a.filter(filt)) == []
    with pytest.raises(ValueError):
        assert [i for i in a if filt(i)] == []


def test_issue_233():
    """
    Make sure hitting a blank line while iterating does not raise IndexError.
    """
    tmp = pybedtools.BedTool._tmp()
    with open(tmp, "w") as fout:
        fout.write(
            dedent(
                """

            chr1\t1\t5

            # chr2\t5\t9
            """
            )
        )
    x = pybedtools.BedTool(tmp)

    # Previously raised IndexError:
    print(x)


def test_issue_246():
    a = pybedtools.BedTool(
        """
        chr1    14831331        14831332        0       name1   A       25      16      0       9       0       0
        chr1    14831623        14831624        0       name2   A       23      16      0       7       0       0
        chr2    7730095 7730096 0       name3   A       20      18      0       2       0       0
        chr2    7735877 7735878 0       name4   A       25      16      0       9       0       0
        """,
        from_string=True,
    )
    b = pybedtools.BedTool(
        """
        chr1    14805135        14882224        geneA   100     +
        """,
        from_string=True,
    )
    ab = a.intersect(b, loj=True)
    assert ab.file_type == "bed"
    print(ab)


def test_issue_251():
    g = pybedtools.example_bedtool("a.bed")
    i = g[0]
    i
    i.fields
    i.attrs

    # previously, this would raise
    # "ValueError: Interval.attrs was not None, but this was a non-GFF Interval
    #
    i.fields


def test_issue_257():
    try:
        import pandas
        import numpy as np
    except ImportError:
        pytest.mark.skip("Pandas not installed; skipping")
    df = pybedtools.example_bedtool("a.bed").to_dataframe()
    df.iloc[-1, -3:] = np.nan
    b = pybedtools.BedTool.from_dataframe(df)
    assert str(b) == fix(
        """
        chr1        1       100     feature1        0.0     +
        chr1        100     200     feature2        0.0     +
        chr1        150     500     feature3        0.0     -
        chr1        900     950     .               .       .
        """
    )


def test_issue_258():
    """
    Non-BED format BedTool objects can still use to_dataframe and use their own
    header
    """
    a = pybedtools.BedTool(
        """
        chr1  1 5
        chr1  5 10
        """,
        from_string=True,
    )
    tmp = pybedtools.BedTool._tmp()
    with open(tmp, "w") as fout:
        fout.write(">chr1\n" "ACACGACTACACTGACTGTGTCGACTAGCACTACGACTGCAGGCATATAC\n")
    b = a.nucleotide_content(fi=tmp)
    df = b.to_dataframe(disable_auto_names=True)
    assert list(df.columns) == [
        "#1_usercol",
        "2_usercol",
        "3_usercol",
        "4_pct_at",
        "5_pct_gc",
        "6_num_A",
        "7_num_C",
        "8_num_G",
        "9_num_T",
        "10_num_N",
        "11_num_oth",
        "12_seq_len",
    ]


def test_issue_303():
    # Issue 303 describes hitting a cap of 253 -b files. Locally I hit a limit
    # at 510, and observe the same on travis-ci.
    #
    # The fix was to check the args in bedtool._wraps, and raise an exception
    # if there's more than 510 filenames provided. Note that it works find with
    # many BedTool objects.

    b = []
    for i in range(1000):
        b.append(
            pybedtools.BedTool(
                "chr1\t{0}\t{1}\tb{0}".format(i, i + 1), from_string=True
            )
        )
    a = pybedtools.example_bedtool("a.bed")

    # Use many BedTool objects; this works
    x = a.intersect(b, wao=True, filenames=True)

    # Try different cutoffs, providing filenames rather than BedTool objects:
    for n in [64, 256, 510]:
        b2 = [i.fn for i in b[:n]]
        try:
            y = a.intersect(b2)

        # If running on a system that supports <510 filenames, we'll get
        # a BEDToolsError, so catch that and report here
        except pybedtools.helpers.BEDToolsError:
            raise ValueError("Hit a limit at {0} files".format(n))

    # Otherwise, too many filenames should raise a pybedtoolsError as detected
    # by the _wraps() function.
    with pytest.raises(pybedtools.helpers.pybedtoolsError):
        y = a.intersect([i.fn for i in b])


def test_issue_291():
    s = "chr1	10	100\n"
    a = pybedtools.BedTool(s, from_string=True)

    # Create a gzipped file identical to a, bypassing the .saveas() mechanism
    tmpgz = pybedtools.BedTool._tmp() + ".gz"
    with gzip.open(tmpgz, "wt") as fout:
        fout.write(s)
    b = pybedtools.BedTool(tmpgz)

    assert a == b

    prefix = pybedtools.BedTool._tmp()

    # save as uncompressed
    c = a.saveas(prefix)

    # extension triggers compressed output from uncompressed input
    d = a.saveas(prefix + ".gz")

    # compressed output from compressed input.
    #
    # Previously this would fail with:
    # UnicodeDecodeError: 'utf-8' codec can't decode byte 0x8b in position 1: invalid start byte
    #
    # The problem was that only the output compression state was being tracked,
    # so compressed input was being opened as uncompressed. Solution was to
    # track input and output compression states separately.
    e = b.saveas(prefix + "x.gz")

    assert a == b == c == d == e


def test_issue_319():
    vrn_file = os.path.join(testdir, "data", "issue319.vcf.gz")
    spliceslop = os.path.join(testdir, "data", "issue319.bed")
    output_bed = os.path.join(testdir, "data", "issue319.out.bed")
    bt = pybedtools.BedTool(vrn_file).intersect(spliceslop, wa=True, header=True, v=True).saveas(output_bed)


def test_issue_333():
    tmp = pybedtools.BedTool._tmp()
    with open(tmp, 'w') as fout:
        pass
    a = pybedtools.BedTool(tmp)

    # Previously would raise EmptyDataError:
    a.to_dataframe()


def test_issue_343():
    def shift_bed(f, shift):
        f.start += shift
        f.stop += shift
        return f

    a = pybedtools.example_bedtool('a.bed')

    # The fix was to ensure that BedTool.remove_invalid() is always working
    # with a file-based BedTool (whcih means calling .saveas() if needed)
    (
        a
        .each(shift_bed, -200)
        .remove_invalid()
        .sort()
    )


def test_issue_345():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = pybedtools.example_bedtool('c.gff')

    z = a.intersect(b=[b.fn,c.fn], C=True, filenames=True)

    # assert " ".join(z._cmds) == f'intersectBed -filenames -b {b.fn} {c.fn} -a {a.fn} -C'
    assert " ".join(z._cmds) == f'intersectBed -a {a.fn} -filenames -b {b.fn} {c.fn} -C'

def test_issue_348():
    i = pybedtools.Interval('chr1', 1, 100, 'feature1', '.', '.', otherfields=['f1'])


def test_issue_355():
    vcf = pybedtools.example_bedtool('v.vcf')
    for line in open(vcf.fn):
        if not line.startswith('#'):
            break
    assert line.split('\t')[1] == '14'
    assert vcf[0].start == 13
