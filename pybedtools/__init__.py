import os
import sys
import subprocess
import tempfile
import urllib2
import copy_reg
import logging
import scripts
from cbedtools import Interval, IntervalFile, overlap, \
    create_interval_from_list, Attributes, \
    MalformedBedLineError, IntervalIterator
import contrib
from _Window import Window
from helpers import get_tempdir, set_tempdir, cleanup, \
    find_tagged, set_bedtools_path, set_samtools_path
import helpers
from bedtool import BedTool
import genome_registry
import stats
from __main__ import main
from version import __version__

import settings

example_files = ['a.bed.', 'b.bed', 'test.fa', 'a.bam']

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.WARNING)
formatter = logging.Formatter('%(name)s [%(levelname)s]: %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


def debug_mode(x):
    """
    Enable debug mode.

    Use debug_mode(True) to show debug log events in the console and to save
    calling info in BedTool objects, and turn it off again with
    debug_mode(False).

    Note that `pybedtools.KEEP_TEMPFILES` will be set as well, so you will need
    to clean up the tempfile directory manually after using debug mode.
    """
    if x:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        _DEBUG = True
        KEEP_TEMPFILES = True
        logger.info(
            'Debug mode enabled.  You may also want to set '
            'pybedtools.KEEP_TEMPFILES=True to prevent automatic deletion '
            'of files upon exit.')
    else:
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)
        _DEBUG = False
        KEEP_TEMPFILES = False
        logger.info('Debug mode disabled')


def check_for_bedtools(program_to_check="intersectBed", force_check=False):
    """
    For backwards compatibility; please use helpers._check_for_bedtools()
    """
    return helpers._check_for_bedtools(program_to_check, force_check)


# Allow Interval objects to be pickled -- required if you want to pass them
# across process boundaries
def interval_constructor(fields):
    return create_interval_from_list(list(fields))


def interval_reducer(interval):
    return interval_constructor, (tuple(interval.fields), )

copy_reg.pickle(Interval, interval_reducer, interval_constructor)


def load_path_config(fn):
    """
    You can use a config file to specify installation paths of various programs
    used by pybedtools.  This can be useful for testing, or using different
    versions of programs.

    `fn` is a config file with the following format.  If an entry is blank,
    then assume it's already on the path. All items must be lowercase::

        [paths]
        bedtools=/tools/BEDTools/bin
        samtools=
        r=
        tabix=

    You only need to specify paths you need to change, so this is a valid file
    that will only specify the path to use for R::

        [paths]
        r=/usr/bin/R-dev

    If `fn` is not a string, then assume it is a dictionary of (program,
    paths). This is used primarily for testing.
    """
    setters = dict(
        bedtools=helpers.set_bedtools_path,
        samtools=helpers.set_samtools_path,
        r=helpers.set_R_path,
        tabix=helpers.set_tabix_path)

    if isinstance(fn, dict):
        for prog, setter in setters.items():
            try:
                path = fn[prog]
                setter(path)
            except KeyError:
                pass

    if isinstance(fn, basestring):
        import ConfigParser
        c = ConfigParser.SafeConfigParser()
        c.read(fn)
        if c.sections() != ['paths']:
            raise ValueError("Invalid path config -- must have "
                             "only one section, [paths].")
        for prog, setter in setters.items():
            try:
                path = c.get('paths', prog)
                setter(path)

            except ConfigParser.NoOptionError:
                pass


def data_dir():
    """
    Returns the data directory that contains example files for tests and
    documentation.
    """
    return os.path.join(os.path.dirname(__file__), 'test', 'data')


def example_filename(fn):
    """
    Return a bed file from the pybedtools examples directory.  Use
    func:`list_example_files` to see a list of files that are included.
    """
    fn = os.path.join(data_dir(), fn)
    if not os.path.exists(fn):
        raise ValueError("%s does not exist" % fn)
    return fn


def example_bedtool(fn):
    """
    Return a bedtool using a bed file from the pybedtools examples directory.
    Use :func:`list_example_files` to see a list of files that are included.
    """
    fn = os.path.join(data_dir(), fn)
    if not os.path.exists(fn):
        raise ValueError("%s does not exist" % fn)
    return BedTool(fn)


def list_example_files():
    """
    Returns a list of files in the examples dir.  Choose one and pass it to
    :func:`example_filename` to get the full path to an example file.

    Example usage:

        >>> choices = list_example_files()
        >>> assert 'a.bed' in choices
        >>> bedfn = example_filename('a.bed')
        >>> mybedtool = BedTool(bedfn)

    """
    candidate_fns = os.listdir(data_dir())
    exts = ('.bed', '.gff', '.gtf', '.bed.gz', '.bam', '.gff.gz')
    valid_fns = [f for f in candidate_fns if f.endswith(exts)]
    return sorted(valid_fns)


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
        ('chr2L', (0, 23011544))
        ('chr2LHet', (0, 368872))
        ('chr2R', (0, 21146708))
        ('chr2RHet', (0, 3288761))
        ('chr3L', (0, 24543557))
        ('chr3LHet', (0, 2555491))
        ('chr3R', (0, 27905053))
        ('chr3RHet', (0, 2517507))
        ('chr4', (0, 1351857))
        ('chrM', (0, 19517))
        ('chrU', (0, 10049037))
        ('chrUextra', (0, 29004656))
        ('chrX', (0, 22422827))
        ('chrXHet', (0, 204112))
        ('chrYHet', (0, 347038))


    """
    try:
        return getattr(genome_registry, genome)
    except AttributeError:
        return get_chromsizes_from_ucsc(genome)


def chromsizes_to_file(chrom_sizes, fn=None):
    """
    Converts a *chromsizes* dictionary to a file.  If *fn* is None, then a
    tempfile is created (which can be deleted with pybedtools.cleanup()).

    Returns the filename.
    """
    if fn is None:
        tmpfn = tempfile.NamedTemporaryFile(prefix='pybedtools.',
                                            suffix='.tmp', delete=False)
        tmpfn = tmpfn.name
        BedTool.TEMPFILES.append(tmpfn)
        fn = tmpfn
    if isinstance(chrom_sizes, basestring):
        chrom_sizes = chromsizes(chrom_sizes)
    fout = open(fn, 'w')
    for chrom, bounds in sorted(chrom_sizes.items()):
        line = chrom + '\t' + str(bounds[1]) + '\n'
        fout.write(line)
    fout.close()
    return fn


def get_chromsizes_from_ucsc(genome, saveas=None, mysql='mysql', timeout=None):
    """
    Download chrom size info for *genome* from UCSC and returns the dictionary.

    If you need the file, then specify a filename with *saveas* (the dictionary
    will still be returned as well).

    If ``mysql`` is not on your path, specify where to find it with
    *mysql=<path to mysql executable>*.

    *timeout* is how long to wait for a response; mostly used for testing.

    Example usage:

        >>> dm3_chromsizes = get_chromsizes_from_ucsc('dm3')
        >>> for i in sorted(dm3_chromsizes.items()):
        ...     print i
        ('chr2L', (0, 23011544))
        ('chr2LHet', (0, 368872))
        ('chr2R', (0, 21146708))
        ('chr2RHet', (0, 3288761))
        ('chr3L', (0, 24543557))
        ('chr3LHet', (0, 2555491))
        ('chr3R', (0, 27905053))
        ('chr3RHet', (0, 2517507))
        ('chr4', (0, 1351857))
        ('chrM', (0, 19517))
        ('chrU', (0, 10049037))
        ('chrUextra', (0, 29004656))
        ('chrX', (0, 22422827))
        ('chrXHet', (0, 204112))
        ('chrYHet', (0, 347038))

    """
    if not internet_on(timeout=timeout):
        raise ValueError('It appears you don\'t have an internet connection '
                         '-- unable to get chromsizes from UCSC')
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
        stdout, stderr = p.communicate()
        if stderr:
            print(stderr)
            print('Commands were:\n')
            print(subprocess.list2cmdline(cmds))

        lines = stdout.splitlines()[1:]
        d = {}
        for line in lines:
            chrom, size = line.split()
            d[chrom] = (0, int(size))

        if saveas is not None:
            chromsizes_to_file(d, saveas)

        return d

    except OSError as err:
        if err.errno == 2:
            raise OSError("Can't find mysql -- if you don't have it "
                          "installed, you'll have to get chromsizes "
                          " manually, or "
                          "specify the path with the 'mysql' kwarg.")
        else:
            raise


def internet_on(timeout=1):
    try:
        response = urllib2.urlopen('http://google.com', timeout=timeout)
        return True
    except urllib2.URLError as err:
        pass
    return False
