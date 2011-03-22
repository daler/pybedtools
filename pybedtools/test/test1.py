import pybedtools
import os, difflib
from nose.tools import *
from pybedtools import bedtool

testdir = os.path.dirname(__file__)

pybedtools.set_tempdir('.')

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
    a = pybedtools.bedtool(os.path.join(testdir, 'a.bed'))
    b = pybedtools.bedtool(os.path.join(testdir, 'b.bed'))
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
    

def decorator_test():
    from pybedtools.bedtool import _returns_bedtool, _help

    @_returns_bedtool()
    def dummy():
        pass
    assert "returns a new bedtool instance" in dummy.__doc__
    
    @_help('intersectBed')
    def dummy2():
        pass
    
    # "-a" ought to be in the help string for intersectBed somewhere....
    assert '-a' in dummy2.__doc__



def test_midpoint():
    a = pybedtools.bedtool("chr1 1 100\nchr5 3000 4000", from_string=True)
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
    assert 'a.bed' in pybedtools.list_example_beds()

    a = pybedtools.example_bed_fn('a.bed')
    assert a == os.path.join(testdir, 'a.bed')
    
    a = pybedtools.example_bedtool('a.bed')
    assert a.fn == os.path.join(testdir, 'a.bed')

    # complain appropriately if nonexistent paths are asked for
    assert_raises(ValueError, pybedtools.example_bed_fn, 'nonexistent')
    assert_raises(ValueError, pybedtools.example_bedtool, 'nonexistent')
    assert_raises(ValueError, pybedtools.set_tempdir, 'nonexistent')

def test_bedtool_creation():
    # make sure we can make a bedtool from a bedtool and that it points to the
    # same file
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.bedtool(a)
    assert b.fn == a.fn
    assert_raises(ValueError, pybedtools.bedtool,'nonexistent.bed')

    # note that *s* has both tabs and spaces....
    s = """
    chr1	1	100	feature1  0	+
    chr1	100	200	feature2  0	+
    chr1	150	500	feature3  0	-
    chr1	900	950	feature4  0	+
    """
    from_string = pybedtools.bedtool(s, from_string=True)

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
    from_string = pybedtools.bedtool(s, from_string=True)
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    
    assert from_string == a
    assert from_string != b


def test_add_subtract():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    assert a.intersect(b,u=True) == (a+b)
    assert a.intersect(b,v=True) == (a-b)

def test_flatten():
    from pybedtools.bedtool import _flatten_list 
    result = _flatten_list([[1,2,3,0,[0,5],9],[100]])
    print result
    assert result == [1, 2, 3, 0, 0, 5, 9, 100]

def test_history_step():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = a.intersect(b)
    d = c.subtract(a)
    
    print d.history
    d.delete_temporary_history(ask=False)
    print d.history

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
    a = pybedtools.bedtool(s, from_string=True)
    b = a.sequence(fi=fi)
    print open(b.seqfn).read()
    assert False, 'Left off here . . . .'
    
    

    
def teardown():
    # always run this!
    pybedtools.cleanup(remove_all=True)
