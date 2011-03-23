import os
import sys
import subprocess
import tempfile

from bedtool import bedtool, get_tempdir, set_tempdir, cleanup, find_tagged


__version__ = '0.2.2dev'

def check_for_bedtools(program_to_check='intersectBed'):
    try:
        p = subprocess.Popen([program_to_check], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError as err:
        if err.errno == 2:
            raise OSError("Please make sure you have installed BEDTools (https://github.com/arq5x/bedtools) and that it's on the path.")

check_for_bedtools()

def data_dir():
    """
    Returns the data directory that contains example files for tests and
    documentation.
    """
    return os.path.join(os.path.dirname(__file__), 'test')

def example_bed_fn(bed):
    """
    Return a bed file from the pybedtools examples directory.  Use
    :func:`list_example_beds` to see a list of files that are included.
    """
    fn = os.path.join(data_dir(), bed)
    if not os.path.exists(fn):
        raise ValueError, "%s does not exist" % fn
    return fn

def example_bedtool(fn):
    """
    Return a bedtool using a bed file from the pybedtools examples directory.
    Use :func:`list_example_beds` to see a list of files that are included.
    """
    fn = os.path.join(data_dir(), fn)
    if not os.path.exists(fn):
        raise ValueError, "%s does not exist" % fn
    return bedtool(fn)

def list_example_beds():
    """
    Returns a list of bed files in the examples dir.  Choose one and pass
    it to :func:`example_bed` to get the full path to an example BED file.

    Example usage::

        >>> choices = list_example_beds()
        >>> bedfn = example_bed(choices[0])
        >>> mybedtool = bedtool(bedfn)


    """
    return sorted([i for i in os.listdir(data_dir()) if os.path.splitext(i)[-1] == '.bed'])


def chromsizes(genome):
    """
    Looks for a *genome* already included in the genome_registry; if not found
    then it looks it up on UCSC.  Returns the dictionary of chromsizes.

    Example usage::
        
        dm3_chromsizes = chromsizes('dm3')

    """
    try:
        return getattr(genome_registry, genome)
    except AttributeError:
        return get_chromsizes_from_ucsc(genome)

def chromsizes_to_file(chromsizes, fn=None):
    """
    Converts a *chromsizes* dictionary to a file.  If *fn* is None, then a
    tempfile is created (which can be deleted with pybedtools.cleanup()).  

    Returns the filename.
    """
    if fn is None:
        tmpfn = tempfile.NamedTemporaryFile(prefix='pybedtools.',suffix='.tmp',delete=False)
        tmpfn = tmpfn.name
        bedtool.TEMPFILES.append(tmpfn)
        fn = tmpfn
    fout = open(fn,'w')
    for chrom, bounds in sorted(chromsizes.items()):
        line = chrom + '\t' + str(bounds[1]) + '\n'
        fout.write(line)
    fout.close()
    return fn

def get_chromsizes_from_ucsc(genome, saveas=None, mysql='mysql'):
    """
    Download chrom size info for *genome* from UCSC and returns the dictionary. 


    If you need the file, then specify a filename with *saveas* (the dictionary
    will still be returned as well).
    
    Example usage::

        a = bedtool('in.bed')
        dm3_chromsizes = a.get_chromsizes_from_ucsc('dm3')

    """
    cmds = [mysql,
            '--user=genome',
            '--host=genome-mysql.cse.ucsc.edu',
            '-A',
            '-e',
            'select chrom, size from %s.chromInfo' % genome]
    try:
        p = subprocess.Popen(cmds,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             bufsize=1)
        stdout,stderr = p.communicate()
        if stderr:
            print stderr
            print 'Commands were:\n'
            print subprocess.list2cmdline(cmds)

        lines = stdout.splitlines()[1:]
        d = {}
        for line in lines:
            chrom,size = line.split()
            d[chrom] = (1,int(size))

        if saveas is not None:
            chromsizes_to_file(d, saveas)

        return d

    except OSError as err:
        if err.errno == 2:
            raise OSError("Can't find mysql -- if you don't have it "
                          "installed, you'll have to get chromsizes manually, or " 
                          "specify the path with the 'mysql' kwarg.")
        else:
            raise
