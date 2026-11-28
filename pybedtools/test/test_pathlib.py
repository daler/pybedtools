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



@pytest.mark.parametrize("make_paths", [list, tuple])
def test_intersect_multiple_pathlib_databases(tmp_path, make_paths):
    import subprocess

    a = pybedtools.example_bedtool("a.bed")
    b = pathlib.Path(pybedtools.example_filename("b.bed"))
    c = tmp_path / "other peaks.bed"
    c.write_text("chr1\t20\t30\tpeak\t0\t+\n")
    expected = subprocess.check_output(
        ["bedtools", "intersect", "-a", a.fn, "-b", str(b), str(c), "-wa", "-wb", "-filenames"]
    )

    result = a.intersect(make_paths((b, c)), wa=True, wb=True, filenames=True)

    assert str(b).encode() in expected
    assert str(c).encode() in expected
    assert pathlib.Path(result.fn).read_bytes() == expected
