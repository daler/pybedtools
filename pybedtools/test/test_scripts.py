import pybedtools
from pybedtools.scripts import annotate
from nose.tools import assert_raises



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
    assert f.fields[-1] == "exon;intron"


def test_annotate_xstream():
    a = pybedtools.example_bedtool('m1.bed')
    b = pybedtools.example_bedtool('mm9.bed12')
    c = annotate.add_xstream(a, b, dist=1000, updown="up")
    assert a.field_count() == c.field_count() - 1
    assert len(a) == len(c)
    d = annotate.add_xstream(c, b, dist=1000, updown="down")
    assert a.field_count() == d.field_count() - 2
