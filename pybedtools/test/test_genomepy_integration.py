from six.moves import builtins
import os
import sys
import pytest


# Make sure all tests import pybedtools + genomepy anew
@pytest.fixture(autouse=True)
def cleanup_imports():
    yield
    del_mods = [mod for mod in sys.modules if "pybedtools" in mod or "genomepy" in mod]
    for mod in del_mods:
        del sys.modules[mod]


def test_genomepy_not_installed():
    import pybedtools

    if "genomepy" in sys.modules:
        del sys.modules["genomepy"]
    genome = "pybedtools/test/data/genome.fa"
    d = pybedtools.helpers.get_chromsizes_from_genomepy(genome)
    assert d is None

    with pytest.raises(OSError):
        pybedtools.chromsizes(genome)
    with pytest.raises(OSError):
        pybedtools.chromsizes("non-existing")


def test_chromsizes_from_genomepy():
    import pybedtools

    if "genomepy" not in sys.modules:
        pytest.skip("genomepy not instlled -- skipping test")

    genome = "pybedtools/test/data/genome.fa"
    try:
        d = pybedtools.helpers.get_chromsizes_from_genomepy(genome)
        print(pybedtools.chromsizes(genome))
        assert d["chr1"] == (0, 10)

        d = pybedtools.chromsizes(genome)
        assert d["chr3"] == (0, 30)

    finally:
        # Make sure all genomepy files get deleted
        fnames = [genome + ext for ext in [".fai", ".sizes"]]
        fnames.append(genome.replace(".fa", ".gaps.bed"))
        for fname in fnames:
            if os.path.exists(fname):
                os.unlink(fname)

    assert None == pybedtools.helpers.get_chromsizes_from_genomepy("non-existing")
