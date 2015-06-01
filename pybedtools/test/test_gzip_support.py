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

def test_gzip():
    # make new gzipped files on the fly
    agz = pybedtools.BedTool._tmp()
    bgz = pybedtools.BedTool._tmp()
    os.system('gzip -c %s > %s' % (pybedtools.example_filename('a.bed'), agz))
    os.system('gzip -c %s > %s' % (pybedtools.example_filename('b.bed'), bgz))
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
