"""
Tests for BedTool.liftover().

Use a tiny hand-written chain file that maps old chr1:0-1000 one-to-one onto
new chrX:100-1100, so every lifted coordinate is simply the old one shifted by
+100 onto chrX.
"""
import os
import shutil
import subprocess

import pytest

import pybedtools
from pybedtools import filenames

needs_liftover = pytest.mark.skipif(
    shutil.which("liftOver") is None, reason="UCSC liftOver not on PATH"
)

# Chain header fields:
#   score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
# followed by ungapped block sizes; the last block line has only a size.
# "t" is the old assembly and "q" is the new one.
CHAIN = """\
chain 1000 chr1 1000 + 0 1000 chrX 2000 + 100 1100 1
1000

"""

# f1 and f2 are in the chain. f3 is on a chromosome the chain doesn't
# cover. f4 is only half covered, so it fails liftOver's default
# -minMatch=0.95.
BED = """\
chr1\t10\t20\tf1\t0\t+
chr1\t500\t600\tf2\t0\t-
chr2\t10\t20\tf3\t0\t+
chr1\t990\t1010\tf4\t0\t+
"""


def teardown_module():
    pybedtools.cleanup()


@pytest.fixture
def chain(tmp_path):
    p = tmp_path / "old_to_new.chain"
    p.write_text(CHAIN)
    return p


@pytest.fixture
def bed():
    return pybedtools.BedTool(BED, from_string=True)


def _records(fn):
    """Non-comment lines of a BED file as lists of fields, used for verifying
    results below."""
    with open(fn) as f:
        return [line.rstrip("\n").split("\t") for line in f if not line.startswith("#")]


@needs_liftover
def test_liftover_maps_and_writes_unmapped(bed, chain, tmp_path):
    # also test spaces in filenames
    d = tmp_path / "dir with space"
    d.mkdir()
    unmapped = d / "unmapped.bed"

    res = bed.liftover(chain, unmapped=unmapped)

    assert _records(res.fn) == [
        ["chrX", "110", "120", "f1", "0", "+"],
        ["chrX", "600", "700", "f2", "0", "-"],
    ]
    assert sorted(r[3] for r in _records(unmapped)) == ["f3", "f4"]


@needs_liftover
def test_liftover_discards_unmapped_by_default(bed, chain):
    res = bed.liftover(chain)
    assert [r[3] for r in _records(res.fn)] == ["f1", "f2"]


@needs_liftover
@pytest.mark.parametrize("args", ["-minMatch=0.4", ["-minMatch=0.4"]])
def test_liftover_args_are_passed(bed, chain, tmp_path, args):
    # f4 is 50% covered, so lowering minMatch should switch it from unmapped to
    # mapped.
    unmapped = tmp_path / "unmapped.bed"
    res = bed.liftover(chain, unmapped=unmapped, liftover_args=args)

    mapped = {r[3]: r for r in _records(res.fn)}
    assert set(mapped) == {"f1", "f2", "f4"}
    assert mapped["f4"][0] == "chrX"
    assert [r[3] for r in _records(unmapped)] == ["f3"]


@needs_liftover
def test_liftover_from_iterator(bed, chain):
    it = pybedtools.BedTool(iter(bed))
    assert not isinstance(it.fn, str)
    res = it.liftover(chain)
    assert [r[3] for r in _records(res.fn)] == ["f1", "f2"]


@needs_liftover
def test_liftover_failure_raises(bed, tmp_path):
    with pytest.raises(subprocess.CalledProcessError) as excinfo:
        bed.liftover(tmp_path / "missing.chain")
    assert excinfo.value.stderr
