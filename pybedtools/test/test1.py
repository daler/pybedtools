import pybedtools
import os, difflib, sys
from nose.tools import assert_raises

testdir = os.path.dirname(__file__)

pybedtools.set_tempdir('.')

def fix(x):
    """
    Replaces spaces with tabs, removes spurious newlines, and lstrip()s each
    line. Makes it really easy to create BED files on the fly for testing and
    checking.
    """
    s = ""
    for i in  x.splitlines():
        i = i.strip()
        if len(i) == 0:
            continue
        i = i.split()
        i = '\t'.join(i)+'\n'
        s += i
    return s

# ----------------------------------------------------------------------------
# Streaming and non-file BedTool tests
# ----------------------------------------------------------------------------
def test_stream():
    """
    Stream and file-based equality, both whole-file and Interval by
    Interval
    """
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = a.intersect(b)

    # make an unwriteable dir...
    orig_tempdir = pybedtools.get_tempdir()
    if os.path.exists('unwriteable'):
        os.system('rm -rf unwriteable')
    os.system('mkdir unwriteable')
    os.system('chmod -w unwriteable')

    # ...set that to the new tempdir
    pybedtools.set_tempdir('unwriteable')

    # this should really not be written anywhere
    d = a.intersect(b, stream=True)

    assert_raises(NotImplementedError, c.__eq__, d)
    d_contents = d.fn.read()
    c_contents = open(c.fn).read()
    assert d_contents == c_contents

    # reconstruct d and check Interval-by-Interval equality
    pybedtools.set_tempdir('unwriteable')
    d = a.intersect(b, stream=True)

    for i,j in zip(c, d):
        assert str(i) == str(j)

    # Now do something similar with GFF files.
    a = pybedtools.example_bedtool('a.bed')
    f = pybedtools.example_bedtool('d.gff')

    # file-based
    pybedtools.set_tempdir(orig_tempdir)
    g1 = f.intersect(a)

    # streaming
    pybedtools.set_tempdir('unwriteable')
    g2 = f.intersect(a, stream=True)

    for i,j in zip(g1, g2):
        assert str(i) == str(j)

    # this was segfaulting at one point, just run to make sure
    g3 = f.intersect(a, stream=True)
    for i in iter(g3):
        print i

    pybedtools.set_tempdir(orig_tempdir)
    os.system('rm -fr unwriteable')

def test_stream_of_stream():
    """
    Second-level streaming using self-intersections
    """
    a = pybedtools.example_bedtool('a.bed')

    # Ensure non-stream and stream equality of self-intersection
    nonstream1 = a.intersect(a, u=True)
    stream1    = a.intersect(a, u=True, stream=True)
    nonstream1_str = str(nonstream1)
    stream1_str    = str(stream1)
    a_str          = str(a)
    assert nonstream1_str == stream1_str == a_str

    # Have to reconstruct stream1 cause it was consumed in the str() call
    nonstream1 = a.intersect(a, u=True)
    stream1    = a.intersect(a, u=True, stream=True)
    nonstream2 = a.intersect(nonstream1, u=True)
    stream2    = a.intersect(stream1, u=True, stream=True)
    nonstream2_str = str(nonstream2)
    stream2_str    = str(stream2)
    assert nonstream2_str == stream2_str == nonstream1_str == stream1_str == a_str

def test_generator():
    """
    Equality of BedTools created from file, iter(), and generator
    """
    # Test creation from file vs 
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.BedTool(iter(a))
    assert str(a) == str(b)

    # Ensure that streams work well too
    b1 = a.intersect(a, stream=True)
    b2 = pybedtools.BedTool((i for i in a)).intersect(a)
    assert str(b1) == str(b2)

#TODO: this hangs -- it's the same as above, but with stream=True on the last
# intersection . . .  is this a Popen issue?
def tst_():
    a = pybedtools.example_bedtool('a.bed')
    b1 = a.intersect(a, stream=True)
    b2 = pybedtools.BedTool((i for i in a)).intersect(a, stream=True)
    assert str(b1) == str(b2)

def test_malformed():
    """
    Malformed BED lines should raise MalformedBedLineError
    """
    a = pybedtools.BedTool("""
    chr1 100 200
    chr1 100 90
    chr1 100 200
    chr1 100 200
    chr1 100 200
    chr1 100 200
    """, from_string=True)
    a_i = iter(a)

    # first feature is OK
    print a_i.next()

    # but next one is not and should raise exception
    assert_raises(pybedtools.MalformedBedLineError, a_i.next)

def test_remove_invalid():
    """
    Remove_invalid() removes invalid lines, track lines, and comments
    """
    a = pybedtools.BedTool("""
    chr1 100 200
    chr1 100 90
    track name='try to break parser'
    chr1 100 200
    chr1 100 200
    chr1 100 200
    #
    chr1 100 200
    """, from_string=True)

    b = a.remove_invalid()

    cleaned = pybedtools.BedTool("""
    chr1 100 200
    chr1 100 200
    chr1 100 200
    chr1 100 200
    chr1 100 200""", from_string=True)

    assert_raises(NotImplementedError, b.__eq__, cleaned)
    assert str(b) == str(cleaned)

def test_create_from_list_long_features():
    """
    Iterator handles extra fields from long features (BED+GFF -wao intersection)
    """
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('c.gff')
    c = a.intersect(b, wao=True, stream=False)
    d = a.intersect(b, wao=True, stream=True)
    for i in d:
        print i

def test_iterator():
    """
    Iterator should ignore non-BED lines
    """
    s = """
    track name="test"


    browser position chrX:1-100
    # comment line
    chrX  1 10
    # more comments
    track name="another"


    """
    a = pybedtools.BedTool(s, from_string=True)
    results = list(a)
    assert str(results[0]) == 'chrX\t1\t10', results

def test_indexing():
    """
    Indexing into BedTools
    """
    a = pybedtools.example_bedtool('a.bed')

    # This is the first line
    interval = pybedtools.Interval('chr1', 1, 100, 'feature1', '0', '+')

    # just to make sure
    assert interval == iter(a).next()

    # test slice behavior
    results = list(a[0:2])
    assert len(results) == 2
    assert results[0] == interval

    # test single-integer indexing
    assert a[0] == interval

    # only slices and integers allowed....
    assert_raises(ValueError, a.__getitem__, 'key')

def test_repr_and_printing():
    """
    Missing files and streams should say so in repr()
    """
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = a+b
    d = a.intersect(b, stream=True)
    os.unlink(c.fn)
    assert 'a.bed' in repr(a)
    assert 'b.bed' in repr(b)
    assert 'MISSING FILE' in repr(c)
    assert 'stream' in repr(d)

def test_file_type():
    """
    Regression test on file_type checks

    Previously file_type was creating a new IntervalFile every time it was
    called; now it's cached so an IntervalFile is only created once per
    BedTool.
    """
    a = pybedtools.example_bedtool('a.bed')
    for i in range(5000):
        a.file_type

# ----------------------------------------------------------------------------
# BEDTools wrapper tests --
#   See test_iter.py, which uses YAML test case definitions, for more complete
#   tests of BEDTools wrapper methods. 
#
#   Here, we assert exception raises and more complicated things that can't be
#   easily described in YAML
# ----------------------------------------------------------------------------

def test_introns():
    a = pybedtools.example_bedtool('mm9.bed12')
    b = pybedtools.BedTool((f for f in a if f.name == "Tcea1,uc007afj.1")).saveas()
    bfeat = iter(b).next()

    bi = b.introns()
    # b[9] is the exonCount from teh bed12 file. there should be
    # b[9] -1 introns assuming no utrs.
    assert len(bi) == int(bfeat[9]) - 1, (len(bi), len(b))

def test_slop():
    """
    Calling slop with no genome should raise ValueError
    """
    a = pybedtools.example_bedtool('a.bed')

    # Make sure it complains if no genome is set
    assert_raises(ValueError, a.slop, **dict(l=100, r=1))

def test_closest():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    r = a.closest(b)
    assert len(r) == len(a)

# TODO: there's enough stuff in here that it's probably worth it to eventually
# make a TestSequenceStuff class
def test_sequence():
    """
    From UCSC:

    chromStart - The starting position of the feature in the chromosome or
    scaffold. The first base in a chromosome is numbered 0.

    chromEnd - The ending position of the feature in the chromosome or
    scaffold. The chromEnd base is not included in the display of the feature.
    For example, the first 100 bases of a chromosome are defined as
    chromStart=0, chromEnd=100, and span the bases numbered 0-99. """

    fi = os.path.join(testdir, 'test.fasta')

    s = """
    chrX 9  16 . . +
    chrX 9  16 . . -
    chrY 1  4  . . +
    chrZ 28 31 . . +
    """

    fasta = """
    >chrX
    AAAAAAAAATGCACTGAAAAAAAAAAAAAAA
    >chrY
    GCTACCCCCCCCCCCCCCCCCCCCCCCCCCC
    >chrZ
    AAAAAAAAAAAAAAAAAAAAAAAAAAAATCT
    """
    a = pybedtools.BedTool(s, from_string=True)
    assert_raises(ValueError, a.save_seqs, ('none',))

    fout = open(fi,'w')
    for line in fasta.splitlines(True):
        fout.write(line.lstrip())
    fout.close()

    # redirect stderr for the call to .sequence(), which reports the creation
    # of an index file
    tmp = open(a._tmp(),'w')
    orig_stderr = sys.stderr
    sys.stderr = tmp

    f = a.sequence(fi=fi)

    sys.stderr = orig_stderr

    assert f.fn == f.fn
    seqs = open(f.seqfn).read()
    print seqs
    expected = """>chrX:9-16
TGCACTG
>chrX:9-16
TGCACTG
>chrY:1-4
CTA
>chrZ:28-31
TCT
"""
    print ''.join(difflib.ndiff(seqs,expected))
    print expected 
    assert seqs == expected

    f = a.sequence(fi=fi,s=True)
    seqs = open(f.seqfn).read()
    expected = """>chrX:9-16(+)
TGCACTG
>chrX:9-16(-)
CAGTGCA
>chrY:1-4(+)
CTA
>chrZ:28-31(+)
TCT
"""
    print seqs
    print expected
    print ''.join(difflib.ndiff(seqs,expected))
    assert seqs == expected

    f = f.save_seqs('deleteme.fa')
    assert open('deleteme.fa').read() == expected
    assert f.print_sequence() == expected
    os.unlink('deleteme.fa')

    fresh_a = pybedtools.BedTool(s, from_string=True)
    assert fresh_a == f

    os.unlink(fi)
    if os.path.exists(fi+'.fai'):
        os.unlink(fi+'.fai')

# ----------------------------------------------------------------------------
# Operator tests
# ----------------------------------------------------------------------------
def test_add_subtract():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    assert a.intersect(b,u=True) == (a+b)
    assert a.intersect(b,v=True) == (a-b)

def test_subset():
    a = pybedtools.example_bedtool('a.bed')
    import random
    random.seed(1)

    s = list(a.random_subset(1).features())
    assert len(s) == 1
    assert isinstance(s[0], pybedtools.Interval)

    s2 = list(a.random_subset(len(a)).features())
    print len(s2)
    assert len(s2) == len(a)

def test_eq():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('a.bed')

    # BedTool to BedTool
    assert a == b

    # BedTool to string
    assert a == """chr1	1	100	feature1	0	+
chr1	100	200	feature2	0	+
chr1	150	500	feature3	0	-
chr1	900	950	feature4	0	+
"""
    # Test not equa on bedtool
    b = pybedtools.example_bedtool('b.bed')
    assert b != a

    # and string
    assert a != "blah"

    # Don't allow testing equality on streams
    c = a.intersect(b, stream=True)
    d = a.intersect(b)
    assert_raises(NotImplementedError, c.__eq__, d)
    assert_raises(NotImplementedError, d.__eq__, c)

    # Test it on iterator, too....
    e = pybedtools.BedTool((i for i in a))
    assert_raises(NotImplementedError, e.__eq__, a)
    assert_raises(NotImplementedError, a.__eq__, e)

    # Make sure that if we force the iterator to be consumed, it is in fact
    # equal
    assert a == str(e)


# ----------------------------------------------------------------------------
# Other BedTool method tests
# ----------------------------------------------------------------------------

def test_count_bed():
    a = pybedtools.example_bedtool('a.bed')
    assert a.count() == 4
    assert len(a) == 4

def test_feature_centers():
    from pybedtools import featurefuncs
    a = pybedtools.BedTool("""
                           chr1 1 100
                           chr5 3000 4000
                           """, from_string=True)
    b = a.each(featurefuncs.center, 1)
    results = list(b.features())

    print results

    assert results[0].start == 50
    assert results[0].stop == 51
    assert results[0].chrom == 'chr1'

    assert results[1].start == 3500
    assert results[1].stop == 3501
    assert results[1].chrom == 'chr5'

def test_bedtool_creation():
    # make sure we can make a bedtool from a bedtool and that it points to the
    # same file
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.BedTool(a)
    assert b.fn == a.fn
    assert_raises(ValueError, pybedtools.BedTool,'nonexistent.bed')

    # note that *s* has both tabs and spaces....
    s = """
    chr1	1	100	feature1  0	+
    chr1	100	200	feature2  0	+
    chr1	150	500	feature3  0	-
    chr1	900	950	feature4  0	+
    """
    from_string = pybedtools.BedTool(s, from_string=True)

    # difflib used here to show a bug where a newline was included when using
    # from_string
    print ''.join(difflib.ndiff(str(from_string), str(a)))

    assert str(from_string) == str(a)

def test_special_methods():
    # note that *s* has both tabs and spaces....
    s = """
    chr1	1	100	feature1  0	+
    chr1	100	200	feature2  0	+
    chr1	150	500	feature3  0	-
    chr1	900	950	feature4  0	+
    """
    from_string = pybedtools.BedTool(s, from_string=True)
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')

    assert from_string == a
    assert from_string != b
    assert not from_string == b
    assert not from_string != a

def test_field_count():
    a = pybedtools.example_bedtool('a.bed')
    assert a.field_count() == 6

def test_repr_and_printing():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = a+b
    os.unlink(c.fn)
    assert 'a.bed' in repr(a)
    assert 'b.bed' in repr(b)
    assert 'MISSING FILE' in repr(c)

    print a.head(1)

def test_cut():
    a = pybedtools.example_bedtool('a.bed')
    c = a.cut([0, 1, 2, 4])
    assert c.field_count() == 4, c

def test_filter():
    a = pybedtools.example_bedtool('a.bed')

    b = a.filter(lambda f: f.length < 100 and f.length > 0)
    assert len(b) == 2

def test_random_intersection():
    # TODO:
    return
    N = 4
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    li = list(a.randomintersection(b, N))
    assert len(li) == N, li

def test_cat():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    b_fn = pybedtools.example_filename('b.bed')
    assert a.cat(b) == a.cat(b_fn)
    expected =  fix("""
    chr1 1   500
    chr1 800 950
    """)
    assert a.cat(b) == expected

    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = a.cat(b, postmerge=False)
    assert len(a) + len(b) == len(c), (len(a), len(b), len(c))


# ----------------------------------------------------------------------------
# Interval tests
# ----------------------------------------------------------------------------

def test_gff_stuff():
    s = """
    chr1  fake  gene 1 100 . + . ID=gene1
    chr1  fake  mRNA 1 100 . + . Name=mRNA1
    chr1  fake  CDS 50 90 . + . other=nothing
    """
    d = pybedtools.BedTool(s, from_string=True)
    f1, f2, f3 = d.features()
    assert f1.name == 'gene1', f1.name
    assert f2.name == 'mRNA1', f2.name
    assert f3.name is None, f3.name

def test_name():
    c = iter(pybedtools.example_bedtool('c.gff')).next()
    assert c.name == "thaliana_1_465_805" , c.name

# ----------------------------------------------------------------------------
# Other non-BedTool tests
# ----------------------------------------------------------------------------

def test_flatten():
    from pybedtools.helpers import _flatten_list 
    result = _flatten_list([[1,2,3,0,[0,5],9],[100]])
    print result
    assert result == [1, 2, 3, 0, 0, 5, 9, 100]

def test_history_step():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = a.intersect(b)
    d = c.subtract(a)

    tag = c.history[0].result_tag
    assert pybedtools.find_tagged(tag) == c

    assert_raises(ValueError, pybedtools.find_tagged, 'nonexistent')


    print d.history
    d.delete_temporary_history(ask=True, raw_input_func=lambda x: 'n')
    assert os.path.exists(a.fn)
    assert os.path.exists(b.fn)
    assert os.path.exists(c.fn)
    assert os.path.exists(d.fn)

    d.delete_temporary_history(ask=True, raw_input_func=lambda x: 'Yes')
    assert os.path.exists(a.fn)
    assert os.path.exists(b.fn)
    assert not os.path.exists(c.fn) # this is the only thing that should change
    assert os.path.exists(d.fn)

    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = a.intersect(b)
    d = c.subtract(a)
    d.delete_temporary_history(ask=False)
    assert os.path.exists(a.fn)
    assert os.path.exists(b.fn)
    assert not os.path.exists(c.fn) # this is the only thing that should change
    assert os.path.exists(d.fn)

def test_kwargs():
    a = pybedtools.example_bedtool('a.bed')
    b = a.intersect(a, s=False)
    c = a.intersect(a)
    assert str(b) == str(c)


# ----------------------------------------------------------------------------
# BAM support tests
# ----------------------------------------------------------------------------
def test_bam_bedtool_creation():
    x = pybedtools.example_bedtool('x.bam')
    a = pybedtools.example_bedtool('a.bed')
    assert x._isbam
    assert not a._isbam

def test_print_abam():
    x = pybedtools.example_bedtool('gdc.bam')
    expected = fix("""
    None	0	chr2L	11	255	5M	*	0	0	CGACA	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	71	255	5M	*	0	0	TTCTC	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	141	255	5M	*	0	0	CACCA	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	151	255	5M	*	0	0	GTTCA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	211	255	5M	*	0	0	AAATA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	71	255	5M	*	0	0	GAGAA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	141	255	5M	*	0	0	TGGTG	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	161	255	5M	*	0	0	GATAA	IIIII	NM:i:0	NH:i:1""")
    print x
    print expected
    assert x == expected

def test_bam_iter():
    x = pybedtools.example_bedtool('gdc.bam')
    s = 'None	0	chr2L	11	255	5M	*	0	0	CGACA	IIIII	NM:i:0	NH:i:1'
    assert str(x[0]) == str(iter(x).next()) == s

def test_bam_stream_bed():
    x = pybedtools.example_bedtool('gdc.bam')
    b = pybedtools.example_bedtool('gdc.gff')
    c = x.intersect(b, u=True, bed=True, stream=True)
    expected = fix("""
    chr2L	70	75	None	255	-
    chr2L	140	145	None	255	-
    chr2L	150	155	None	255	-
    chr2L	210	215	None	255	+
    chr2L	70	75	None	255	+
    chr2L	140	145	None	255	+
    chr2L	160	165	None	255	+
    """)
    assert c == expected

def test_bam_stream_bam():
    x = pybedtools.example_bedtool('gdc.bam')
    b = pybedtools.example_bedtool('gdc.gff')
    c = x.intersect(b, u=True, stream=True)
    expected = fix("""
    None	16	chr2L	71	255	5M	*	0	0	TTCTC	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	141	255	5M	*	0	0	CACCA	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	151	255	5M	*	0	0	GTTCA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	211	255	5M	*	0	0	AAATA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	71	255	5M	*	0	0	GAGAA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	141	255	5M	*	0	0	TGGTG	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	161	255	5M	*	0	0	GATAA	IIIII	NM:i:0	NH:i:1""")
    assert str(c) == expected

# TODO: stream-of-stream doesn't work yet for BAMs
def test_bam_stream_bam_stream():
    x = pybedtools.example_bedtool('gdc.bam')
    b = pybedtools.example_bedtool('gdc.gff')
    c = x.intersect(b, u=True, stream=True)
    expected = fix("""
    None	16	chr2L	71	255	5M	*	0	0	TTCTC	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	141	255	5M	*	0	0	CACCA	IIIII	NM:i:0	NH:i:1
    None	16	chr2L	151	255	5M	*	0	0	GTTCA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	211	255	5M	*	0	0	AAATA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	71	255	5M	*	0	0	GAGAA	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	141	255	5M	*	0	0	TGGTG	IIIII	NM:i:0	NH:i:1
    None	0	chr2L	161	255	5M	*	0	0	GATAA	IIIII	NM:i:0	NH:i:1""")   
    d = c.intersect(b)
    print d
    assert str(d) == expected

def test_bam_interval():
    x = pybedtools.example_bedtool('x.bam')
    assert x[0].chrom == 'chr2L'
    assert x[0].start == 9329L
    assert x[0][3] == '9330'
    assert x[0].stop == 9365L
    assert len(x[0][9]) == len(x[0]) == 36

def test_bam_regression():
    # Regression test:  with extra fields, the first item in x.bam was being
    # parsed as gff (cause not ==13 fields).  This does a check to prevent that
    # from happening again.
    x = pybedtools.example_bedtool('x.bam')
    assert x[0].file_type == 'sam'
    assert x[0].chrom == 'chr2L'

def teardown():
    # always run this!
    pybedtools.cleanup(remove_all=True)
