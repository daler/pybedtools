import os
import sys
import subprocess
import tempfile


from helpers import get_tempdir, set_tempdir, cleanup, find_tagged
from bedtool import BedTool
import genome_registry


__version__ = '0.2.2dev'

# Registry of files in the test dir that should be accessible to users via example_bedtool
example_files = ['a.bed.','b.bed','test.fa', 'a.bam']

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

def example_filename(fn):
    """
    Return a bed file from the pybedtools examples directory.  Use
    :func:`list_example_files` to see a list of files that are included.
    """
    fn = os.path.join(data_dir(), fn)
    if not os.path.exists(fn):
        raise ValueError, "%s does not exist" % fn
    return fn

def example_bedtool(fn):
    """
    Return a bedtool using a bed file from the pybedtools examples directory.
    Use :func:`list_example_files` to see a list of files that are included.
    """
    fn = os.path.join(data_dir(), fn)
    if not os.path.exists(fn):
        raise ValueError, "%s does not exist" % fn
    return BedTool(fn)

def list_example_files():
    """
    Returns a list of files in the examples dir.  Choose one and pass it to
    :func:`example_file_fnl` to get the full path to an examplefile.

    Example usage:

        >>> choices = list_example_files()
        >>> assert 'a.bed' in choices
        >>> bedfn = example_filename('a.bed')
        >>> mybedtool = BedTool(bedfn)
        >>> print mybedtool
        chr1 1   100 feature1 0 +
        chr1 100 200 feature2 0 +
        chr1 150 500 feature3 0 -
        chr1 900 950 feature4 0 +
        <BLANKLINE>

    """
    return sorted([i for i in os.listdir(data_dir()) if os.path.splitext(i)[-1] == '.bed'])

def chromsizes(genome):
    """
    Looks for a *genome* already included in the genome registry; if not found
    then it looks it up on UCSC.  Returns the dictionary of chromsize tuples
    where each tuple has (start,stop).

    Chromsizes are described as (start, stop) tuples to allow randomization
    within specified regions; e. g., you can make a chromsizes dictionary that
    represents the extent of a tiling array.

    Example usage:

        >>> dm3_chromsizes = chromsizes('dm3')
        >>> for i in sorted(dm3_chromsizes.items()):
        ...     print i
        ('chr2L', (1, 23011544))
        ('chr2LHet', (1, 368872))
        ('chr2R', (1, 21146708))
        ('chr2RHet', (1, 3288761))
        ('chr3L', (1, 24543557))
        ('chr3LHet', (1, 2555491))
        ('chr3R', (1, 27905053))
        ('chr3RHet', (1, 2517507))
        ('chr4', (1, 1351857))
        ('chrM', (1, 19517))
        ('chrU', (1, 10049037))
        ('chrUextra', (1, 29004656))
        ('chrX', (1, 22422827))
        ('chrXHet', (1, 204112))
        ('chrYHet', (1, 347038))


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
        BedTool.TEMPFILES.append(tmpfn)
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

    Example usage:

        >>> dm3_chromsizes = get_chromsizes_from_ucsc('dm3')
        >>> for i in sorted(dm3_chromsizes.items()):
        ...     print i
        ('chr2L', (1, 23011544))
        ('chr2LHet', (1, 368872))
        ('chr2R', (1, 21146708))
        ('chr2RHet', (1, 3288761))
        ('chr3L', (1, 24543557))
        ('chr3LHet', (1, 2555491))
        ('chr3R', (1, 27905053))
        ('chr3RHet', (1, 2517507))
        ('chr4', (1, 1351857))
        ('chrM', (1, 19517))
        ('chrU', (1, 10049037))
        ('chrUextra', (1, 29004656))
        ('chrX', (1, 22422827))
        ('chrXHet', (1, 204112))
        ('chrYHet', (1, 347038))

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

if __name__ == "__main__":
    print 'Running tests...'
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
