from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pybedtools.test.tfuncs as tfuncs

import os
import pybedtools
import gzip

setup = tfuncs.setup
teardown = tfuncs.teardown

def _make_temporary_gzip(bed_filename):
    gz_filename = pybedtools.BedTool._tmp()
    with gzip.open(gz_filename, 'w') as out_:
        with open(bed_filename, 'r') as in_:
            out_.writelines(in_)
    return gz_filename

def test_gzip():
    # make new gzipped files on the fly
    agz = _make_temporary_gzip(pybedtools.example_filename('a.bed'))
    bgz = _make_temporary_gzip(pybedtools.example_filename('b.bed'))

    agz = pybedtools.BedTool(agz)
    bgz = pybedtools.BedTool(bgz)
    assert agz.file_type == bgz.file_type == 'bed'
    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    assert a.intersect(b) == agz.intersect(bgz) == a.intersect(bgz) == agz.intersect(b)


def test_gzipped_output():
    expected = pybedtools.BedTool._tmp()
    fn = pybedtools.example_filename('a.bed')
    os.system('gzip -c {fn} > {expected}'.format(**locals()))
    obs = pybedtools.example_bedtool('a.bed').saveas(compressed=True)
    assert gzip.open(obs.fn).read() == gzip.open(expected).read()
