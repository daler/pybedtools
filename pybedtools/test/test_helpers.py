import pybedtools
import six
import sys
import os, difflib
from .tfuncs import setup_module, teardown_module, testdir, test_tempdir, unwriteable
import pytest


def fix(x):
    """
    Replaces spaces with tabs, removes spurious newlines, and lstrip()s each
    line. Makes it really easy to create BED files on the fly for testing and
    checking.
    """
    s = ""
    for i in x.splitlines():
        i = i.strip()
        if len(i) == 0:
            continue
        i = i.split()
        i = "\t".join(i) + "\n"
        s += i
    return s


def test_isBAM():
    bam = pybedtools.example_filename("x.bam")
    notabam = pybedtools.example_filename("a.bed")
    open("tiny.txt", "w").close()
    assert pybedtools.helpers.isBAM(bam)
    assert not pybedtools.helpers.isBAM(notabam)
    assert not pybedtools.helpers.isBAM("tiny.txt")
    os.unlink("tiny.txt")


def test_cleanup():
    """
    make sure the tempdir and cleanup work
    """
    # assert os.path.abspath(pybedtools.get_tempdir()) == os.path.abspath('.')

    # make a fake tempfile, not created during this pybedtools session
    pybedtools.cleanup()

    testfn = os.path.join(test_tempdir, "pybedtools.TESTING.tmp")
    os.system("touch %s" % testfn)
    assert os.path.exists(testfn)

    # make some temp files
    a = pybedtools.BedTool(os.path.join(testdir, "data", "a.bed"))
    b = pybedtools.BedTool(os.path.join(testdir, "data", "b.bed"))
    c = a.intersect(b)

    # after standard cleanup, c's fn should be gone but the fake one still
    # there...
    pybedtools.cleanup(verbose=True)
    assert os.path.exists(testfn)
    assert not os.path.exists(c.fn)

    # Unless we force the removal of all temp files.
    pybedtools.cleanup(remove_all=True)
    assert not os.path.exists(testfn)

    # a.fn and b.fn better be there still!
    assert os.path.exists(a.fn)
    assert os.path.exists(b.fn)


# Test disabled, since we are no longer checking for pre-v2.15 versions by using
# separate program calls
#
# TODO: remove this once v2.15 support is fully deprecated.
#
# def test_bedtools_check():
#     # this should run fine (especially since we've already imported pybedtools)
#     pybedtools.check_for_bedtools()
#
#     assert_raises(OSError, pybedtools.check_for_bedtools, **dict(program_to_check='nonexistent', force_check=True))


def test_version_check():
    pybedtools.helpers._check_for_bedtools(
        force_check=True, override="bedtools v2.28", verbose=True
    )
    assert pybedtools.settings.bedtools_version == [2, 28]
    pybedtools.helpers._check_for_bedtools(
        force_check=True, override="bedtools debian/2.28.0+dfsg-2-dirty", verbose=True
    )
    assert pybedtools.settings.bedtools_version == [2, 28, 0]


def test_call():
    tmp = os.path.join(pybedtools.get_tempdir(), "test.output")
    from pybedtools.helpers import call_bedtools, BEDToolsError

    with pytest.raises(BEDToolsError):
        call_bedtools(*(["intersectBe"], tmp))

    a = pybedtools.example_bedtool("a.bed")

    # momentarily redirect stderr to file so the error message doesn't spew all
    # over the place when testing
    orig_stderr = sys.stderr
    sys.stderr = open(a._tmp(), "w")
    sys.stderr = orig_stderr

    pybedtools.set_bedtools_path("nonexistent")
    a = pybedtools.example_bedtool("a.bed")
    with pytest.raises(NotImplementedError):
        a.intersect(a)
    pybedtools.set_bedtools_path()
    a = pybedtools.example_bedtool("a.bed")
    assert a.intersect(a, u=True) == a


def test_chromsizes():
    with pytest.raises(OSError):
        pybedtools.get_chromsizes_from_ucsc(
            "dm3", mysql="wrong path", fetchchromsizes="wrongtoo"
        )
    with pytest.raises(ValueError):
        pybedtools.get_chromsizes_from_ucsc("dm3", timeout=0)
    try:

        print(pybedtools.chromsizes("dm3"))
        print(pybedtools.get_chromsizes_from_ucsc("dm3"))
        assert pybedtools.chromsizes("dm3") == pybedtools.get_chromsizes_from_ucsc(
            "dm3"
        )

        hg17 = pybedtools.chromsizes("hg17")

        assert hg17["chr1"] == (0, 245522847)

        fn = pybedtools.chromsizes_to_file(hg17, fn="hg17.genome")
        expected = "chr1\t245522847\n"
        results = open(fn).readline()
        print(results)
        assert expected == results

        # make sure the tempfile version works, too
        fn = pybedtools.chromsizes_to_file(hg17, fn=None)
        expected = "chr1\t245522847\n"
        results = open(fn).readline()
        print(results)
        assert expected == results

        with pytest.raises(OSError):
            pybedtools.get_chromsizes_from_ucsc(
                **dict(genome="hg17", mysql="nonexistent", fetchchromsizes="missing")
            )

        os.unlink("hg17.genome")
    except OSError:
        sys.stdout.write("mysql error -- test for chromsizes from UCSC didn't run")


def test_ff_center():
    from pybedtools.featurefuncs import center

    a = pybedtools.example_bedtool("a.bed")
    b = a.each(center, width=10)
    expected = fix(
        """
    chr1	45	55	feature1	0	+
    chr1	145	155	feature2	0	+
    chr1	320	330	feature3	0	-
    chr1	920	930	feature4	0	+"""
    )
    assert str(b) == expected


def test_getting_example_beds():
    assert "a.bed" in pybedtools.list_example_files()

    a_fn = pybedtools.example_filename("a.bed")
    assert a_fn == os.path.join(testdir, "data", "a.bed")

    a = pybedtools.example_bedtool("a.bed")
    assert a.fn == os.path.join(testdir, "data", "a.bed")

    # complain appropriately if nonexistent paths are asked for
    e = FileNotFoundError if six.PY3 else ValueError
    with pytest.raises(e):
        pybedtools.example_filename("nonexistent")
    with pytest.raises(e):
        pybedtools.example_bedtool("nonexistent")
    with pytest.raises(e):
        pybedtools.set_tempdir("nonexistent")


def teardown():
    # always run this!
    pybedtools.cleanup(remove_all=True)
