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
    scale_fn = bam + ".scale"
    if os.path.exists(scale_fn) and not force:
        for line in open(scale_fn):
            if line.startswith("#"):
                continue
            readcount = float(line.strip())
            return readcount

    cmds = ["samtools", "view", "-c", "-F", "0x4", bam]
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode:
        raise ValueError("samtools says: %s" % stderr)

    readcount = float(stdout)

    # write to file so the next time you need the lib size you can access
    # it quickly
    if not os.path.exists(scale_fn):
        fout = open(scale_fn, "w")
        fout.write(str(readcount) + "\n")
        fout.close()
    return readcount


def bedgraph_to_bigwig(bedgraph, genome, output):
    genome_file = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))
    cmds = ["bedGraphToBigWig", bedgraph.fn, genome_file, output]
    try:
        p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
    except FileNotFoundError:
        raise FileNotFoundError(
            "bedGraphToBigWig was not found on the path. This is an external "
            "tool from UCSC which can be downloaded from "
            "http://hgdownload.soe.ucsc.edu/admin/exe/. Alternatatively, use "
            "`conda install ucsc-bedgraphtobigwig`"
        )

    if p.returncode:
        raise ValueError(
            "cmds: %s\nstderr:%s\nstdout:%s" % (" ".join(cmds), stderr, stdout)
        )
    return output


def bigwig_to_bedgraph(fn, chrom=None, start=None, end=None, udcDir=None):
    cmds = ["bigWigToBedGraph", fn]
    if chrom is not None:
        cmds.extend(["-chrom", chrom])
    if start is not None:
        cmds.extend(["-start", start])
    if end is not None:
        cmds.extend(["-end", end])
    if udcDir is not None:
        cmds.extend(["-udcDir", udcDir])

    outfn = pybedtools.BedTool._tmp()
    cmds.append(outfn)

    try:
        p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
    except FileNotFoundError:
        raise FileNotFoundError(
            "bigWigToBedGraph was not found on the path. This is an external "
            "tool from UCSC which can be downloaded from "
            "http://hgdownload.soe.ucsc.edu/admin/exe/. Alternatatively, use "
            "`conda install ucsc-bedgraphtobigwig`"
        )
    if p.returncode:
        raise ValueError(
            "cmds: %s\nstderr:%s\nstdout:%s" % (" ".join(cmds), stderr, stdout)
        )
    return pybedtools.BedTool(outfn)


def wig_to_bigwig(wig, genome, output):
    genome_file = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))
    cmds = ["wigToBigWig", wig.fn, genome_file, output]

    try:
        p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
    except FileNotFoundError:
        raise FileNotFoundError(
            "bigWigToBedGraph was not found on the path. This is an external "
            "tool from UCSC which can be downloaded from "
            "http://hgdownload.soe.ucsc.edu/admin/exe/. Alternatatively, use "
            "`conda install ucsc-bedgraphtobigwig`"
        )
    if p.returncode:
        raise ValueError(
            "cmds: %s\nstderr:%s\nstdout:%s" % (" ".join(cmds), stderr, stdout)
        )
    return output


def bam_to_bigwig(bam, genome, output, scale=False):
    """
    Given a BAM file `bam` and assembly `genome`, create a bigWig file scaled
    such that the values represent scaled reads -- that is, reads per million
    mapped reads.

    (Disable this scaling step with scale=False; in this case values will
    indicate number of reads)

    Assumes that `bedGraphToBigWig` from UCSC tools is installed; see
    http://genome.ucsc.edu/goldenPath/help/bigWig.html for more details on the
    format.
    """
    genome_file = pybedtools.chromsizes_to_file(pybedtools.chromsizes(genome))
    kwargs = dict(bg=True, split=True, g=genome_file)
    if scale:
        readcount = mapped_read_count(bam)
        _scale = 1 / (readcount / 1e6)
        kwargs["scale"] = _scale
    x = pybedtools.BedTool(bam).genome_coverage(**kwargs)
    cmds = ["bedGraphToBigWig", x.fn, genome_file, output]
    try:
        p = subprocess.Popen(
            cmds,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        stdout, stderr = p.communicate()
    except FileNotFoundError:
        raise FileNotFoundError(
            "bedGraphToBigWig was not found on the path. This is an external "
            "tool from UCSC which can be downloaded from "
            "http://hgdownload.soe.ucsc.edu/admin/exe/. Alternatatively, use "
            "`conda install ucsc-bedgraphtobigwig`"
        )

    if p.returncode and "bedSort" in stderr:
        print("BAM header was not sorted; sorting bedGraph")
        y = x.sort()
        cmds[1] = y.fn
        try:
            p = subprocess.Popen(
                cmds,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
            stdout, stderr = p.communicate()
        except FileNotFoundError:
            raise FileNotFoundError(
                "bedSort was not found on the path. This is an external "
                "tool from UCSC which can be downloaded from "
                "http://hgdownload.soe.ucsc.edu/admin/exe/. Alternatatively, use "
                "`conda install ucsc-bedgraphtobigwig`"
            )

    if p.returncode:
        raise ValueError(
            "cmds: %s\nstderr: %s\nstdout: %s" % (" ".join(cmds), stderr, stdout)
        )
