"""
when bugs are identified post-release, put tests here to make sure they don't
happen again
"""

import pybedtools
import pybedtools.featurefuncs


def test_midpoint():
    """
    regression test for #98
    """
    a = """chr1 3052874 3053149
    chr1 3333690 3333915
    chr1 3472838 3473382
    chr1 3639053 3639356
    """

    def nothing(f):
        return f

    input_bed = pybedtools.BedTool(
        a, from_string=True).saveas("test_input.bed")

    for func in [pybedtools.featurefuncs.midpoint, pybedtools.featurefuncs.center, nothing]:
        input_bed_mid = input_bed.each(func)
        assert len(input_bed_mid) == 4
