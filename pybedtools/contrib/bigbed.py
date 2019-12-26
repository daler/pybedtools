import os
import subprocess
import six
import pybedtools


def bigbed(
    x,
    genome,
    output,
    blockSize=256,
    itemsPerSlot=512,
    bedtype=None,
    _as=None,
    unc=False,
    tab=False,
):
    """
    Converts a BedTool object to a bigBed format and returns the new filename.

    `x` is a BedTool object

    `genome` is an assembly string

    `output` is the name of the bigBed file to create.

    Other args are passed to bedToBigBed.  In particular, `bedtype` (which
    becomes the "-type=" argument) is automatically handled for you if it is
    kept as the default None.

    Assumes that a recent version of bedToBigBed from UCSC is on the path.
    """
    if isinstance(x, six.string_types):
        x = pybedtools.BedTool(x)
    if not isinstance(x.fn, six.string_types):
        x = x.saveas()
    chromsizes = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))
    if bedtype is None:
        bedtype = "bed%s" % x.field_count()
    cmds = [
        "bedToBigBed",
        x.fn,
        chromsizes,
        output,
        "-blockSize=%s" % blockSize,
        "-itemsPerSlot=%s" % itemsPerSlot,
        "-type=%s" % bedtype,
    ]
    if unc:
        cmds.append("-unc")
    if tab:
        cmds.append("-tab")
    if _as:
        cmds.append("-as=%s" % _as)
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode:
        raise ValueError(
            "cmds: %s\nstderr:%s\nstdout:%s" % (" ".join(cmds), stderr, stdout)
        )

    return output


def bigbed_to_bed(fn, chrom=None, start=None, end=None, maxItems=None):
    cmds = ["bigBedToBed", fn]
    if chrom is not None:
        cmds.extend(["-chrom", chrom])
    if start is not None:
        cmds.extend(["-start", start])
    if end is not None:
        cmds.extend(["-end", end])
    if maxItems is not None:
        cmds.extend(["-maxItems", maxItems])

    outfn = pybedtools.BedTool._tmp()
    cmds.append(outfn)

    p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode:
        raise ValueError(
            "cmds: %s\nstderr:%s\nstdout:%s" % (" ".join(cmds), stderr, stdout)
        )
    return pybedtools.BedTool(outfn)
