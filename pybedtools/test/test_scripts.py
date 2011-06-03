import pybedtools
from pybedtools.scripts import annotate, venn_mpl, venn_gchart
from nose.tools import assert_raises
import os


def test_annotate_main():
    # exits after printing help when sys.argv is not as it should be.
    assert_raises(SystemExit, annotate.main)

def test_annotate_closest():
    a = pybedtools.example_bedtool('m1.bed')
    b = pybedtools.example_bedtool('mm9.bed12')
    c = annotate.add_closest(a, b)
    assert len(a) == len(c), (len(a), len(c), str(c))
    assert a.field_count() == c.field_count() - 2
    # in this test-case, the final column should be exon;intron
    # since m1 completely contains both an exon and an intron.
    f = iter(c).next()
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
    here = os.path.dirname(__file__)
    expected = open(os.path.join(here, 'mpl-expected.png')).read()

    pybedtools.bedtool.random.seed(1)
    a = pybedtools.example_bedtool('rmsk.hg18.chr21.small.bed')
    b = a.random_subset(100).shuffle(genome='hg19', seed=1)
    b = b.cat(a.random_subset(100))
    c = a.random_subset(200).shuffle(genome='hg19', seed=2)
    c = c.cat(b.random_subset(100))

    venn_mpl.venn_mpl(a=a.fn, b=b.fn, c=c.fn, colors=['r','b','g'], outfn='out.png', labels=['a','b','c'])
    assert open(outfn).read() == expected

    os.unlink(outfn)



def test_venn_gchart():
    here = os.path.dirname(__file__)
    expected = open(os.path.join(here, 'gchart-expected.png')).read()


    pybedtools.bedtool.random.seed(1)
    a = pybedtools.example_bedtool('rmsk.hg18.chr21.small.bed')
    b = a.random_subset(100).shuffle(genome='hg19', seed=1)
    b = b.cat(a.random_subset(100))
    c = a.random_subset(200).shuffle(genome='hg19', seed=2)
    c = c.cat(b.random_subset(100))
    colors='00FF00,FF0000,0000FF'
    outfn = 'out.png'
    labels = 'a,b,c'

    expected_data = {'chco': '00FF00,FF0000,0000FF', 'chd': 't:1.0,0.197,0.297,0.102,0.058,0.1,0.058', 'chs': '300x300', 'cht': 'v', 'chdl': 'a|b|c'}

    data = venn_gchart.venn_gchart(a=a.fn,
                            b=b.fn,
                            c=c.fn,
                            colors=colors.split(','),
                            labels=labels.split(','),
                            size='300x300')

    print data
    assert data == expected_data

    venn_gchart.gchart(data, outfn)

    assert open(outfn).read() == expected
    os.unlink(outfn)

def test_venn_mpl_main():
    assert_raises(SystemExit, venn_mpl.main)

def test_venn_gchart_main():
    assert_raises(SystemExit, venn_gchart.main)

def teardown():
    pybedtools.cleanup()
