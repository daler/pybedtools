import os
import pathlib

import six

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
    if six.PY2:
        with pytest.raises(ValueError):
            pybedtools.BedTool(path)
    else:
        with pytest.raises(FileNotFoundError):
            pybedtools.BedTool(path)
