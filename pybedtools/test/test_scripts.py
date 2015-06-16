import pybedtools
from .tfuncs import setup, teardown
from pybedtools.scripts import annotate, venn_mpl, venn_gchart
from nose.tools import assert_raises, assert_equal
from nose.plugins.attrib import attr
from nose.plugins.skip import SkipTest
import os
import sys


def test_annotate_main():
    # exits after printing help when sys.argv is not as it should be.
    orig_stderr = sys.stderr
    sys.stderr = open('annotmp','w')
    assert_raises(SystemExit, annotate.main)
    sys.stderr = orig_stderr
    os.unlink('annotmp')

def test_annotate_closest():
    a = pybedtools.example_bedtool('m1.bed')
    b = pybedtools.example_bedtool('mm9.bed12')
    c = annotate.add_closest(a, b)
    assert len(a) == len(c), (len(a), len(c), str(c))
    assert a.field_count() == c.field_count() - 2
    # in this test-case, the final column should be exon;intron
    # since m1 completely contains both an exon and an intron.
    f = next(iter(c))
    # waiting for fix to bedtools:
    #assert f.fields[-1] == "exon;intron", f.fields[-1]


def test_annotate_xstream():
    a = pybedtools.example_bedtool('m1.bed')
    b = pybedtools.example_bedtool('mm9.bed12')
    c = annotate.add_xstream(a, b, dist=1000, updown="up")
    assert a.field_count() == c.field_count() - 1
    assert len(a) == len(c)
    d = annotate.add_xstream(c, b, dist=1000, updown="down")
    assert a.field_count() == d.field_count() - 2


def test_venn_mpl():
    """
    compares output image to expected
    """
    try:
        import matplotlib
    except ImportError:
        import sys
        sys.stderr.write('Need matplotlib installed to test venn_mpl')
        return

    here = os.path.dirname(__file__)
    expected_fn = os.path.join(here, 'mpl-expected.png')

    original = pybedtools.example_bedtool('rmsk.hg18.chr21.small.bed').sort().merge()
    a = pybedtools.BedTool(original[:300]).saveas()
    b = pybedtools.BedTool(original[:20]).saveas().cat(pybedtools.BedTool(original[400:500]).saveas())
    c = pybedtools.BedTool(original[15:30]).saveas().cat(pybedtools.BedTool(original[450:650]).saveas())

    outfn = 'mplout.png'
    venn_mpl.venn_mpl(a=a.fn, b=b.fn, c=c.fn, colors=['r','b','g'], outfn=outfn, labels=['a','b','c'])

    # On a different machine, the created image is not visibly different but is
    # numerically different.  Not sure what a reasonable tolerance is, but this
    # seems to work for now....
    o = matplotlib.image.imread(outfn)
    e = matplotlib.image.imread(expected_fn)

    TOLERANCE = 200
    SUM = abs((o - e).sum())
    assert SUM < TOLERANCE, SUM

    os.unlink(outfn)

def test_venn_gchart_data_is_correct():
    original = pybedtools.example_bedtool('rmsk.hg18.chr21.small.bed').sort().merge()
    a = pybedtools.BedTool(original[:300]).saveas()
    b = pybedtools.BedTool(original[:20]).saveas().cat(pybedtools.BedTool(original[400:500]).saveas())
    c = pybedtools.BedTool(original[15:30]).saveas().cat(pybedtools.BedTool(original[450:650]).saveas())

    colors='00FF00,FF0000,0000FF'
    labels = 'a,b,c'

    expected_data = {'chco': '00FF00,FF0000,0000FF',
                     'chd': 't:1.0,0.4,0.7167,0.0667,0.05,0.1833,0.0167',
                     'chs': '300x300',
                     'cht': 'v',
                     'chdl': 'a|b|c'}

    data = venn_gchart.venn_gchart(a=a.fn,
                            b=b.fn,
                            c=c.fn,
                            colors=colors.split(','),
                            labels=labels.split(','),
                            size='300x300')

    for key in expected_data.keys():
        e = expected_data[key]
        o = data[key]
        assert_equal(e, o, 'Key:{!r}\nExpected:{!r}\nObserved:{!r}\n'.format(key, e, o))


def test_venn_gchart_png_is_saved_correctly():
    from nose.plugins.skip import SkipTest
    raise SkipTest('Small differences between fonts in PNG are sometimes produced and need to be accounted for')
    here = os.path.dirname(__file__)
    expected = open(os.path.join(here, 'gchart-expected.png')).read()

    original = pybedtools.example_bedtool('rmsk.hg18.chr21.small.bed').sort().merge()
    a = pybedtools.BedTool(original[:300]).saveas()
    b = pybedtools.BedTool(original[:20]).saveas().cat(pybedtools.BedTool(original[400:500]).saveas())
    c = pybedtools.BedTool(original[15:30]).saveas().cat(pybedtools.BedTool(original[450:650]).saveas())

    colors='00FF00,FF0000,0000FF'
    outfn = 'gchart_out.png'
    labels = 'a,b,c'

    data = venn_gchart.venn_gchart(a=a.fn,
                            b=b.fn,
                            c=c.fn,
                            colors=colors.split(','),
                            labels=labels.split(','),
                            size='300x300')

    venn_gchart.gchart(data, outfn)
    assert_equal(open(outfn).read(), expected)
    os.unlink(outfn)

def test_venn_mpl_main():
    orig_stderr = sys.stderr
    sys.stderr = open('mpltmp','w')
    assert_raises(SystemExit, venn_mpl.main)
    sys.stderr = orig_stderr
    os.unlink('mpltmp')

def test_venn_gchart_main():
    orig_stderr = sys.stderr
    sys.stderr = open('gcharttmp','w')
    assert_raises(SystemExit, venn_gchart.main)
    sys.stderr = orig_stderr
    os.unlink('gcharttmp')
