import pybedtools
import os, difflib
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


def test_cleanup():
    """
    make sure the tempdir and cleanup work
    """
    assert os.path.abspath(pybedtools.get_tempdir()) == os.path.abspath('.')

    # make a fake tempfile, not created during this pybedtools session
    testfn = 'pybedtools.TESTING.tmp'
    os.system('touch %s' % testfn)
    assert os.path.exists(testfn)

    # make some temp files
    a = pybedtools.BedTool(os.path.join(testdir, 'data', 'a.bed'))
    b = pybedtools.BedTool(os.path.join(testdir, 'data', 'b.bed'))
    c = a.intersect(b)

    # after standard cleanup, c's fn should be gone but the fake one still
    # there...
    pybedtools.cleanup(verbose=True)
    assert os.path.exists(testfn)
    assert not os.path.exists(c.fn)

    # Unless we force the removal of all temp files.
    pybedtools.cleanup(remove_all=True)
    assert not os.path.exists(testfn)

    # a.fn and b.fn better be there still!
    assert os.path.exists(a.fn)
    assert os.path.exists(b.fn)

def test_decorators():
    from pybedtools.helpers import _returns_bedtool, _help

    @_returns_bedtool()
    def dummy():
        pass
    assert "returns a new bedtool instance" in dummy.__doc__

    @_help('intersectBed')
    def dummy2():
        pass

    # "-a" ought to be in the help string for intersectBed somewhere....
    assert '-a' in dummy2.__doc__

def test_bedtools_check():
    # this should run fine (especially since we've already imported pybedtools)
    pybedtools.check_for_bedtools()

    # but this should crap out
    assert_raises(OSError, pybedtools.check_for_bedtools, **dict(program_to_check='nonexistent',))

def test_call():
    tmp = os.path.join(pybedtools.get_tempdir(), 'test.output')
    from pybedtools.helpers import call_bedtools, BEDToolsError
    assert_raises(BEDToolsError, call_bedtools, *(['intersectBe'], tmp))

"""
# TODO: test for connection + mysql
def test_chromsizes():
    assert pybedtools.chromsizes('dm3') == pybedtools.get_chromsizes_from_ucsc('dm3')
    
    hg17 = pybedtools.chromsizes('hg17')

    assert hg17['chr1'] == (1,245522847)
   
    fn = pybedtools.chromsizes_to_file(hg17, fn='hg17.genome')
    expected = 'chr1\t245522847\n'
    results = open(fn).readline()
    print results
    assert expected == results

    # make sure the tempfile version works, too
    fn = pybedtools.chromsizes_to_file(hg17, fn=None)
    expected = 'chr1\t245522847\n'
    results = open(fn).readline()
    print results
    assert expected == results

    assert_raises(OSError, pybedtools.get_chromsizes_from_ucsc, **dict(genome='hg17', mysql='nonexistent'))
   
    os.unlink('hg17.genome')
"""

def test_ff_center():
    from pybedtools.featurefuncs import center
    a = pybedtools.example_bedtool('a.bed')
    b = a.each(center, width=10)
    expected = fix("""
    chr1	45	55	feature1	0	+
    chr1	145	155	feature2	0	+
    chr1	320	330	feature3	0	-
    chr1	920	930	feature4	0	+""")
    assert str(b) == expected



def test_stream():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = a.intersect(b)
    d = a.intersect(b, stream=True)
    assert_raises(NotImplementedError, c.__eq__,d)
    d_contents = d.fn.read()
    c_contents = open(c.fn).read()
    assert d_contents == c_contents

    c = a.intersect(b)
    d = a.intersect(b, stream=True)

    for i,j in zip(c, d):
        assert str(i) == str(j)

def test_stream_gen():
    # these should run
    a = pybedtools.example_bedtool('a.bed')
    f = pybedtools.example_bedtool('d.gff')
    g1 = f.intersect(a)
    g2 = f.intersect(a, stream=True)
    for i,j in zip(g1, g2):
        assert str(i) == str(j)

def test_stream_of_stream():
    a = pybedtools.example_bedtool('a.bed')
    stream1  = a.intersect(a, stream=True)
    stream2 = a.intersect(stream1, stream=True)
    s1 = str(stream1)
    s2 = str(stream2)
    print 'stream1'
    print s1
    print 'stream2'
    print s2
    assert s1 == s2

def test_generator():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.BedTool(iter(a))
    expected = str(a)
    observed = str(b)
    print expected
    print observed
    assert expected == observed

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

def test_count_bed():
    a = pybedtools.example_bedtool('a.bed')
    assert a.count() == 4
    assert len(a) == 4

def test_feature_centers():
    return # TODO: fix this
    a = pybedtools.BedTool("""
                           chr1 1 100
                           chr5 3000 4000
                           """, from_string=True)
    b = a.feature_centers(1)
    results = list(b.features())

    print results
    
    assert results[0].start == 50
    assert results[0].stop == 51
    assert results[0].chr == 'chr1'
    
    assert results[1].start == 3500
    assert results[1].stop == 3501
    assert results[1].chr == 'chr5'

def test_getting_example_beds():
    assert 'a.bed' in pybedtools.list_example_files()

    a_fn = pybedtools.example_filename('a.bed')
    assert a_fn == os.path.join(testdir, 'data', 'a.bed')
    
    a = pybedtools.example_bedtool('a.bed')
    assert a.fn == os.path.join(testdir, 'data', 'a.bed')

    # complain appropriately if nonexistent paths are asked for
    assert_raises(ValueError, pybedtools.example_filename, 'nonexistent')
    assert_raises(ValueError, pybedtools.example_bedtool, 'nonexistent')
    assert_raises(ValueError, pybedtools.set_tempdir, 'nonexistent')

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

def test_add_subtract():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    assert a.intersect(b,u=True) == (a+b)
    assert a.intersect(b,v=True) == (a-b)

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
    
    fout = open(fi,'w')
    for line in fasta.splitlines(True):
        fout.write(line.lstrip())
    fout.close()

    b = a.sequence(fi=fi)
    assert b.fn == a.fn
    seqs = open(b.seqfn).read()
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
    
    b = a.sequence(fi=fi,s=True)
    seqs = open(b.seqfn).read()
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
       
    os.unlink(fi)
    if os.path.exists(fi+'.fai'):
        os.unlink(fi+'.fai')

def test_iterator():
    # makes sure we're ignoring non-feature lines
    
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

def test_repr_and_printing():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = a+b
    os.unlink(c.fn)
    assert 'a.bed' in repr(a)
    assert 'b.bed' in repr(b)
    assert 'MISSING FILE' in repr(c)

    print a.head(1)

def test_intersect():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    assert a.intersect(b.fn) == a.intersect(b)


    # straight-up
    expected = fix("""
    chr1 155 200 feature2 0 +
    chr1 155 200 feature3 0 -
    chr1 900 901 feature4 0 +
    """)
    assert str(a.intersect(b)) == expected
    
    # a that have b
    expected = fix("""
    chr1 100 200 feature2 0 +
    chr1 150 500 feature3 0 -
    chr1 900 950 feature4 0 +
    """)
    assert str(a.intersect(b,u=True)) == expected
    
    # stranded straight-up
    expected = fix("""
    chr1 155 200 feature3 0 -
    chr1 900 901 feature4 0 +
    """)
    assert str(a.intersect(b,s=True)) == expected

    # stranded a that have b
    expected = fix("""
    chr1 150 500 feature3 0 -
    chr1 900 950 feature4 0 +
    """)
    assert str(a.intersect(b, u=True, s=True)) == expected

    # a with no b
    expected = fix("""
    chr1 1 100 feature1 0 +
    """)
    assert str(a.intersect(b, v=True)) == expected

    # stranded a with no b
    expected = fix("""
    chr1 1   100 feature1 0 +
    chr1 100 200 feature2 0 +
    """)
    assert str(a.intersect(b, v=True, s=True)) == expected



    
    

def test_subtract():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')

    # plain 'old subtract
    results = str(a.subtract(b))
    expected = fix("""
    chr1	1	100	feature1	0	+
    chr1	100	155	feature2	0	+
    chr1	150	155	feature3	0	-
    chr1	200	500	feature3	0	-
    chr1	901	950	feature4	0	+""")
    assert results == expected

    # strand-specific
    results = str(a.subtract(b,s=True))
    print results
    expected = fix("""
    chr1	1	100	feature1	0	+
    chr1	100	200	feature2	0	+
    chr1	150	155	feature3	0	-
    chr1	200	500	feature3	0	-
    chr1	901	950	feature4	0	+""")
    assert results == expected

    # the difference in f=0.2 and f=0.1 is in feature5 of b.  Since feature2
    # and feature3 of a overlap, it's seeing the 'a' feature as a 399-bp
    # feature (chr1:100-500), and feature5 of 'b' overlaps this by 44 bp.
    #
    # So the threshold fraction should be
    #
    #   44/399 = 0.1103 
    #
    # f > 0.1103 should return no subtractions, because nothing in b overlaps by that much.
    # However, 0.12 doesn't work; need to go to 0.13 . . .
    results = str(a.subtract(b,s=True,f=0.13))
    expected = fix("""
    chr1	1	100	feature1	0	+
    chr1	100	200	feature2	0	+
    chr1	150	500	feature3	0	-
    chr1	900	950	feature4	0	+""")
    assert results == expected
    
    # f < 0.1103, so should get a subtraction
    results = str(a.subtract(b,s=True,f=0.1))
    print results
    expected = fix("""
    chr1	1	100	feature1	0	+
    chr1	100	200	feature2	0	+
    chr1	150	155	feature3	0	-
    chr1	200	500	feature3	0	-
    chr1	900	950	feature4	0	+""")
    assert results == expected

def test_slop():
    a = pybedtools.example_bedtool('a.bed')

    results = str(a.slop(g=pybedtools.chromsizes('hg19'), b=100))
    expected = fix("""
    chr1	0	200	feature1	0	+
    chr1	0	300	feature2	0	+
    chr1	50	600	feature3	0	-
    chr1	800	1050	feature4	0	+
    """)
    assert results == expected

    results = str(a.slop(g=pybedtools.chromsizes('hg19'), l=100, r=1))
    expected = fix("""
    chr1	0	101	feature1	0	+
    chr1	0	201	feature2	0	+
    chr1	50	501	feature3	0	-
    chr1	800	951	feature4	0	+
    """)
    assert results == expected
    
    
    # Make sure it complains if no genome is set
    assert_raises(ValueError, a.slop, **dict(l=100, r=1))

    # OK, so set one...
    a.set_chromsizes(pybedtools.chromsizes('hg19'))

    # Results should be the same as before.
    results = str(a.slop(l=100, r=1))
    expected = fix("""
    chr1	0	101	feature1	0	+
    chr1	0	201	feature2	0	+
    chr1	50	501	feature3	0	-
    chr1	800	951	feature4	0	+
    """)
    print results
    assert results == expected

def test_merge():
    a = pybedtools.example_bedtool('a.bed')
    results = str(a.merge())
    expected = fix("""
    chr1	1	500
    chr1	900	950
    """)
    assert results == expected

    results = str(a.merge(s=True))
    expected = fix("""
    chr1	1	200	+
    chr1	900	950	+
    chr1	150	500	-
    """)
    assert results == expected

    b = pybedtools.example_bedtool('b.bed')
    results = str(b.merge(d=700))
    expected = fix("""
    chr1 155 901 
    """)
    print results
    assert results == expected

def test_closest():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    r = a.closest(b)
    assert len(r) == len(a)

def test_cat():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = a.cat(b, postmerge=False)
    assert len(a) + len(b) == len(c), (len(a), len(b), len(c))

def test_field_count():
    a = pybedtools.example_bedtool('a.bed')
    assert a.field_count() == 6

def test_cut():
    a = pybedtools.example_bedtool('a.bed')
    c = a.cut([0, 1, 2, 4])
    assert c.field_count() == 4, c

def test_name():
    c = iter(pybedtools.example_bedtool('c.gff')).next()
    assert c.name == "thaliana_1_465_805" , c.name

def test_filter():
    a = pybedtools.example_bedtool('a.bed')

    b = a.filter(lambda f: f.length < 100 and f.length > 0)
    assert len(b) == 2

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

def test_random_intersection():
    # TODO:
    return
    N = 4
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    li = list(a.randomintersection(b, N))
    assert len(li) == N, li

def teardown():
    # always run this!
    pybedtools.cleanup(remove_all=True)


