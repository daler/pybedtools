from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import os
from nose.tools import assert_equal, assert_list_equal
import tempfile
import pybedtools.test.tfuncs as tfuncs

import pybedtools
import gzip

setup = tfuncs.setup
teardown = tfuncs.teardown

def _make_temporary_gzip(bed_filename):
    """
    Make a gzip file on the fly
    :param bed_filename: Filename of bed file to gzip
    :return: filename of gzipped file
    """
    orig_suffix = pybedtools.settings.tempfile_suffix
    pybedtools.settings.tempfile_suffix = '.gz'
    gz_filename = pybedtools.BedTool._tmp()
    pybedtools.settings.tempfile_suffix = orig_suffix
    with gzip.open(gz_filename, 'wb') as out_:
        with open(bed_filename, 'rb') as in_:
            out_.writelines(in_)
    return gz_filename

def test_gzipped_file_types_are_bed():
    agz = _make_temporary_gzip(pybedtools.example_filename('a.bed'))

    agz = pybedtools.BedTool(agz)
    assert_equal('bed', agz.file_type)

def test_gzipped_files_can_be_intersected():
    agz = _make_temporary_gzip(pybedtools.example_filename('a.bed'))
    bgz = _make_temporary_gzip(pybedtools.example_filename('b.bed'))

    agz = pybedtools.BedTool(agz)
    bgz = pybedtools.BedTool(bgz)

    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.example_bedtool('b.bed')
    assert a.intersect(b) == agz.intersect(bgz) == a.intersect(bgz) == agz.intersect(b)

def test_gzipped_files_are_iterable_as_normal():
    agz = _make_temporary_gzip(pybedtools.example_filename('a.bed'))
    agz = pybedtools.BedTool(agz)
    a = pybedtools.example_bedtool('a.bed')
    for i in agz:
        print(i)
    assert_list_equal(list(a), list(agz))

def test_str_representation_of_gzipped_files_is_the_same_as_normal():
    agz = _make_temporary_gzip(pybedtools.example_filename('a.bed'))
    agz = pybedtools.BedTool(agz)
    a = pybedtools.example_bedtool('a.bed')
    assert_equal(str(a), str(agz))

def test_head_of_gzipped_files_is_the_same_as_normal():
    agz = _make_temporary_gzip(pybedtools.example_filename('a.bed'))
    agz = pybedtools.BedTool(agz)
    a = pybedtools.example_bedtool('a.bed')
    assert_equal(agz.head(), a.head())

def test_gzipped_output():
    _filename = pybedtools.example_filename('a.bed')
    compressed_file = pybedtools.BedTool(_filename).saveas(compressed=True)

    # Open gzipped file in text mode
    with gzip.open(compressed_file.fn, 'rt') as gf:
        uncompressed_content = gf.read()

    with open(_filename) as f:
        original_content = f.read()

    assert_equal(original_content, uncompressed_content)

def test_gzipping_is_default_when_extension_is_dot_gz():
    _filename = pybedtools.example_filename('a.bed')
    with open(_filename) as f:
        expected_content = f.read()

    __, temp_filename = tempfile.mkstemp(suffix='.gz')
    try:
        bedtool = pybedtools.BedTool(_filename)
        bedtool.saveas(fn=temp_filename)

        with gzip.open(temp_filename, 'rt') as gf:
            # gzip will fail next line if file is not gzipped
            actual_content = gf.read()

        assert_equal(expected_content, actual_content)
    finally:
        if os.path.isfile(temp_filename):
            os.unlink(temp_filename)

def test_gzipping_can_be_turned_off_even_for_dot_gz():
    _filename = pybedtools.example_filename('a.bed')
    with open(_filename) as f:
        expected_content = f.read()

    __, temp_filename = tempfile.mkstemp(suffix='.gz')
    try:
        bedtool = pybedtools.BedTool(_filename)
        bedtool.saveas(fn=temp_filename, compressed=False)

        with open(temp_filename) as non_gz_f:
            # actual content will be jumbled if non_gz_f is unset
            actual_content = non_gz_f.read()

        assert_equal(expected_content, actual_content)
    finally:
        if os.path.isfile(temp_filename):
            os.unlink(temp_filename)

