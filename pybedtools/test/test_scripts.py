import pybedtools
from nose.tools import assert_raises



def test_annotate_main():
    # exits after printing help when sys.argv is not as it should be.
    assert_raises(SystemExit, pybedtools.scripts.annotate.main)

def test_annotate_closest():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = pybedtools.scripts.annotate.add_closest(a, b)
    assert len(a) == len(c), (len(a), len(c), str(c))
    assert a.field_count() == c.field_count() - 2


def test_annotate_xstream():
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    c = pybedtools.scripts.annotate.add_xstream(a, b, dist=1000, updown="up")
    assert a.field_count() == c.field_count() - 1
    assert len(a) == len(c)
    d = pybedtools.scripts.annotate.add_xstream(c, b, dist=1000, updown="down")
    assert a.field_count() == d.field_count() - 2
