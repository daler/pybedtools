import os
import pathlib

import pybedtools
import pytest


def test_pathlib_base():
    file = "a.bed"
    fn = os.path.join(pybedtools.filenames.data_dir(), file)
    path = pathlib.PurePath(fn)
    assert pybedtools.BedTool(path).fn == fn


def test_pathlib_derived():
    file = "a.bed"
    fn = os.path.join(pybedtools.filenames.data_dir(), file)
    path = pathlib.Path(fn)
    assert pybedtools.BedTool(path).fn == fn


# this may be unnecessary, as the check is performed after str conversion
def test_pathlib_nonexistent_file():
    fn = os.path.join(pybedtools.filenames.data_dir(), "this_file_is_missing")
    path = pathlib.Path(fn)
    with pytest.raises(FileNotFoundError):
        pybedtools.BedTool(path)


def test_pathlib_on_methods():
    # Above functions test the creation of BedTool from Path objects; here we
    # test arbitrary methods to test that Path objs are converted to str when
    # composing the relevant command.
    fn_a = os.path.join(pybedtools.filenames.data_dir(), "a.bed")
    fn_b = os.path.join(pybedtools.filenames.data_dir(), "b.bed")
    p_a = pathlib.Path(fn_a)
    p_b = pathlib.Path(fn_b)
    bt = pybedtools.BedTool(fn_a)
    res = bt.intersect(p_b)

