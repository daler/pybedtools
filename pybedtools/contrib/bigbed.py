import os
import subprocess
import pybedtools


def bigbed(x, genome, output, blockSize=256, itemsPerSlot=512, bedFields=None, _as=None, unc=False):
    if isinstance(x, basestring):
        x = pybedtools.BedTool(x)
    if not isinstance(x.fn, basestring):
        x = x.saveas()
    chromsizes = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))
    if bedFields is None:
        bedFields = x.field_count()
    cmds = [
        'bedToBigBed',
        x.fn,
        chromsizes,
        output,
        '-blockSize=%s' % blockSize,
        '-itemsPerSlot=%s' % itemsPerSlot,
        '-bedFields=%s' % bedFields
    ]
    if unc:
        cmds.append('-unc')

    os.system(' '.join(cmds))
    return output

