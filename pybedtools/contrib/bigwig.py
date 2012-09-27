"""
Module to help create scaled bigWig files from BAM
"""
import pybedtools
import os
import subprocess


def mapped_read_count(bam, force=False):
    """
    Scale is cached in a bam.scale file containing the number of mapped reads.
    Use force=True to override caching.
    """
    scale_fn = bam + '.scale'
    if os.path.exists(scale_fn) and not force:
        for line in open(scale_fn):
            if line.startswith('#'):
                continue
            readcount = float(line.strip())
            return readcount

    cmds = ['samtools',
            'view',
            '-c',
            '-F', '0x4',
            bam]
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if stderr:
        raise ValueError('samtools says: %s' % stderr)

    readcount = float(stdout)

    # write to file so the next time you need the lib size you can access
    # it quickly
    if not os.path.exists(scale_fn):
        fout = open(scale_fn, 'w')
        fout.write(str(readcount) + '\n')
        fout.close()
    return readcount


def bam_to_bigwig(bam, genome, output):
    genome_file = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))
    readcount = mapped_read_count(bam)
    scale = 1 / (readcount / 1e6)
    x = pybedtools.BedTool(bam)\
        .genome_coverage(bg=True, scale=scale, split=True, g=genome_file)
    cmds = [
        'bedGraphToBigWig',
        x.fn,
        genome_file,
        output]
    os.system(' '.join(cmds))
