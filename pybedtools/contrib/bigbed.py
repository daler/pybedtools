import os
import subprocess
import pybedtools


def bigbed(x, genome, output, blockSize=256, itemsPerSlot=512, bedtype=None, _as=None, unc=False, tab=False):
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
    if isinstance(x, basestring):
        x = pybedtools.BedTool(x)
    if not isinstance(x.fn, basestring):
        x = x.saveas()
    chromsizes = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))
    if bedtype is None:
        bedtype = 'bed%s' % x.field_count()
    cmds = [
        'bedToBigBed',
        x.fn,
        chromsizes,
        output,
        '-as=%s' % _as,
        '-blockSize=%s' % blockSize,
        '-itemsPerSlot=%s' % itemsPerSlot,
        '-type=%s' % bedtype
    ]
    if unc:
        cmds.append('-unc')
    if tab:
        cmds.append('-tab')
    p = subprocess.Popen(cmds)
    stdout, stderr = p.communicate()
    if p.returncode:
        raise ValueError

    return output
