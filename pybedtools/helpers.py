from __future__ import print_function
import sys
import os
import tempfile
import subprocess
import random
import string
import glob
import struct
import atexit
import six
from six.moves import urllib
from . import cbedtools
from . import settings
from . import filenames
from . import genome_registry
from .logger import logger
from .cbedtools import create_interval_from_list as _create_interval_from_list


BUFSIZE = 1

_tags = {}

def create_interval_from_list(x):
    y = [i.encode('UTF-8') for i in x]
    for i in y:
        assert type(i) == bytes
    return _create_interval_from_list(y)

def _check_for_bedtools(program_to_check='intersectBed', force_check=False):
    """
    Checks installation as well as version (based on whether or not "bedtools
    intersect" works, or just "intersectBed")
    """
    if settings._bedtools_installed and not force_check:
        return True

    try:
        p = subprocess.Popen(
            [os.path.join(settings._bedtools_path, 'bedtools'),
             settings._prog_names[program_to_check]],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        settings._bedtools_installed = True
        settings._v_2_15_plus = True

    except (OSError, KeyError) as err:

        try:
            p = subprocess.Popen(
                [os.path.join(settings._bedtools_path, program_to_check)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            settings._bedtools_installed = True
            settings._v_2_15_plus = False

        except OSError as err:
            if err.errno == 2:
                if settings._bedtools_path:
                    add_msg = "(tried path '%s')" % settings._bedtools_path
                else:
                    add_msg = ""
                raise OSError("Please make sure you have installed BEDTools"
                              "(https://github.com/arq5x/bedtools) and that "
                              "it's on the path. %s" % add_msg)


def _check_for_tabix():
    try:
        p = subprocess.Popen(
            [os.path.join(settings._tabix_path, 'tabix')],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        settings._tabix_installed = True
    except OSError:
        if settings._tabix_path:
            add_msg = "(tried path '%s')" % settings._tabix_path
        else:
            add_msg = ""
        raise ValueError(
            'Please install tabix and ensure it is on your path %s'
            % add_msg)


def _check_for_samtools():
    try:
        p = subprocess.Popen(
            [os.path.join(settings._samtools_path, 'samtools')],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        settings._samtools_installed = True
    except OSError:
        if settings._samtools_path:
            add_msg = "(tried path '%s')" % settings._samtools_path
        else:
            add_msg = ""
        raise ValueError(
            'Please install samtools and ensure it is on your path %s'
            % add_msg)


def _check_for_R():
    try:
        p = subprocess.Popen(
            [os.path.join(settings._R_path, 'R'), '--version'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        settings._R_installed = True
    except OSError:
        if settings._R_path:
            add_msg = "(tried path '%s')" % settings._R_path
        else:
            add_msg = ""
        raise ValueError(
            'Please install R and ensure it is on your path %s' % add_msg)


class Error(Exception):
    """Base class for this module's exceptions"""
    pass


class BEDToolsError(Error):
    def __init__(self, cmd, msg):
        self.cmd = str(cmd)
        self.msg = str(msg)

    def __str__(self):
        m = '\nCommand was:\n\n\t' + self.cmd + '\n' + \
            '\nError message was:\n' + self.msg
        return m


def isBGZIP(fn):
    """
    Reads a filename to see if it's a BGZIPed file or not.
    """
    header_str = open(fn).read(15)
    if len(header_str) < 15:
        return False

    header = struct.unpack_from('BBBBiBBHBBB', header_str)

    id1, id2, cm, flg, mtime, xfl, os_, xlen, si1, si2, slen = header
    if (id1 == 31) and (id2 == 139) and (cm == 8) and (flg == 4) and \
       (si1 == 66) and (si2 == 67) and (slen == 2):
        return True
    return False


def isBAM(fn):
    if not isBGZIP(fn):
        return False

    # Need to differentiate between BAM and plain 'ol BGZIP. Try reading header
    # . . .
    if not settings._samtools_installed:
        _check_for_samtools()

    cmds = ['samtools', 'view', '-H', fn]
    try:

        # Silence the output, we want to check the return code
        with open(os.devnull, "w") as out:
            subprocess.check_call(cmds, stdout=out, stderr=out)
        return True

    except subprocess.CalledProcessError:
        # Non-0 return code, it means we have an error
        return False

    except OSError:
        raise OSError(
            'SAMtools (http://samtools.sourceforge.net/) '
            'needs to be installed for BAM support')


def find_tagged(tag):
    """
    Returns the bedtool object with tagged with *tag*.  Useful for tracking
    down bedtools you made previously.
    """
    for key, item in _tags.items():
        try:
            if item._tag == tag:
                return item
        except AttributeError:
            pass
    raise ValueError('tag "%s" not found' % tag)


def _flatten_list(x):
    nested = True
    while nested:
        check_again = False
        flattened = []

        for element in x:
            if isinstance(element, list):
                flattened.extend(element)
                check_again = True
            else:
                flattened.append(element)
        nested = check_again
        x = flattened[:]
    return x


def set_tempdir(tempdir):
    """
    Set the directory for temp files.

    Useful for clusters that use a /scratch partition rather than a /tmp dir.
    Convenience function to simply set tempfile.tempdir.
    """
    if not os.path.exists(tempdir):
        errstr = 'The tempdir you specified, %s, does not exist' % tempdir
        raise ValueError(errstr)
    tempfile.tempdir = tempdir


def get_tempdir():
    """
    Gets the current tempdir for the module.
    """
    return tempfile.gettempdir()


def cleanup(verbose=False, remove_all=False):
    """
    Deletes all temp files from the current session (or optionally *all* \
            sessions)

    If *verbose*, reports what it's doing

    If *remove_all*, then ALL files matching "pybedtools.*.tmp" in the temp dir
    will be deleted.
    """
    if settings.KEEP_TEMPFILES:
        return
    for fn in filenames.TEMPFILES:
        if verbose:
            print('removing', fn)
        if os.path.exists(fn):
            os.unlink(fn)
    if remove_all:
        fns = glob.glob(os.path.join(get_tempdir(), 'pybedtools.*.tmp'))
        for fn in fns:
            os.unlink(fn)


def _version_2_15_plus_names(prog_name):
    if not settings._bedtools_installed:
        _check_for_bedtools()
    if not settings._v_2_15_plus:
        return [prog_name]
    try:
        prog_name = settings._prog_names[prog_name]
    except KeyError:
        if prog_name in settings._new_names:
            pass
        raise BEDToolsError(
            prog_name, prog_name + 'not a recognized BEDTools program')
    return [os.path.join(settings._bedtools_path, 'bedtools'), prog_name]


def call_bedtools(cmds, tmpfn=None, stdin=None, check_stderr=None):
    """
    Use subprocess.Popen to call BEDTools and catch any errors.

    Output goes to *tmpfn*, or, if None, output stays in subprocess.PIPE and
    can be iterated over.

    *stdin* is an optional file-like object that will be sent to
    subprocess.Popen.

    Prints some useful help upon getting common errors.

    *check_stderr* is a function that takes the stderr string as input and
    returns True if it's OK (that is, it's not really an error).  This is
    needed, e.g., for calling fastaFromBed which will report that it has to
    make a .fai for a fasta file.
    """
    input_is_stream = stdin is not None
    output_is_stream = tmpfn is None

    _orig_cmds = cmds[:]
    cmds = []
    cmds.extend(_version_2_15_plus_names(_orig_cmds[0]))
    cmds.extend(_orig_cmds[1:])

    try:
        # coming from an iterator, sending as iterator
        if input_is_stream and output_is_stream:
            logger.debug(
                'helpers.call_bedtools(): input is stream, output is '
                'stream')
            logger.debug(
                'helpers.call_bedtools(): cmds=%s', ' '.join(cmds))
            p = subprocess.Popen(cmds,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 stdin=subprocess.PIPE,
                                 bufsize=BUFSIZE)
            for line in stdin:
                p.stdin.write(line)

            # This is important to prevent deadlocks
            p.stdin.close()

            output = p.stdout
            stderr = None

        # coming from an iterator, writing to file
        if input_is_stream and not output_is_stream:
            logger.debug(
                'helpers.call_bedtools(): input is stream, output is file')
            logger.debug(
                'helpers.call_bedtools(): cmds=%s', ' '.join(cmds))
            outfile = open(tmpfn, 'w')
            p = subprocess.Popen(cmds,
                                 stdout=outfile,
                                 stderr=subprocess.PIPE,
                                 stdin=subprocess.PIPE,
                                 bufsize=BUFSIZE)
            if hasattr(stdin, 'read'):
                stdout, stderr = p.communicate(stdin.read())
            else:
                for item in stdin:
                    p.stdin.write(item)
                stdout, stderr = p.communicate()
            output = tmpfn
            outfile.close()

        # coming from a file, sending as iterator
        if not input_is_stream and output_is_stream:
            logger.debug(
                'helpers.call_bedtools(): input is filename, '
                'output is stream')
            logger.debug(
                'helpers.call_bedtools(): cmds=%s', ' '.join(cmds))
            p = subprocess.Popen(cmds,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 bufsize=BUFSIZE)
            output = p.stdout
            stderr = None

        # file-to-file
        if not input_is_stream and not output_is_stream:
            logger.debug(
                'helpers.call_bedtools(): input is filename, output '
                'is filename (%s)', tmpfn)
            logger.debug(
                'helpers.call_bedtools(): cmds=%s', ' '.join(cmds))
            outfile = open(tmpfn, 'w')
            p = subprocess.Popen(cmds,
                                 stdout=outfile,
                                 stderr=subprocess.PIPE,
                                 bufsize=BUFSIZE)
            stdout, stderr = p.communicate()
            output = tmpfn
            outfile.close()

        # Check if it's OK using a provided function to check stderr. If it's
        # OK, dump it to sys.stderr so it's printed, and reset it to None so we
        # don't raise an exception
        if check_stderr is not None:
            if check_stderr(stderr):
                sys.stderr.write(stderr)
                stderr = None

        if stderr:
            raise BEDToolsError(subprocess.list2cmdline(cmds), stderr)

    except (OSError, IOError) as err:
        print('%s: %s' % (type(err), os.strerror(err.errno)))
        print('The command was:\n\n\t%s\n' % subprocess.list2cmdline(cmds))

        problems = {
            2: ('* Did you spell the command correctly?',
                '* Do you have BEDTools installed and on the path?'),
            13: ('* Do you have permission to write '
                 'to the output file ("%s")?' % tmpfn,),
            24: ('* Too many files open -- please submit '
                 'a bug report so that this can be fixed',)
        }

        print('Things to check:')
        print('\n\t' + '\n\t'.join(problems[err.errno]))
        raise OSError('See above for commands that gave the error')

    return output


def set_bedtools_path(path=""):
    """
    Explicitly set path to `BEDTools` installation dir.

    If BEDTools is not available on your system path, specify the path to the
    dir containing the BEDTools executables (intersectBed, subtractBed, etc)
    with this function.

    To reset and use the default system path, call this function with no
    arguments or use path="".
    """
    settings._bedtools_path = path


def set_samtools_path(path=""):
    """
    Explicitly set path to `samtools` installation dir.

    If samtools is not available on the path, then it can be explicitly
    specified here.

    Use path="" to reset to default system path.
    """
    settings._samtools_path = path


def set_tabix_path(path=""):
    """
    Explicitly set path to `tabix` installation dir.

    If tabix is not available on the path, then it can be explicitly
    specified here.

    Use path="" to reset to default system path.
    """
    settings._tabix_path = path


def set_R_path(path=""):
    """
    Explicitly set path to `R` installation dir.

    If R is not available on the path, then it can be explicitly
    specified here.

    Use path="" to reset to default system path.
    """
    settings._R_path = path


def _check_sequence_stderr(x):
    """
    If stderr created by fastaFromBed starst with 'index file', then don't
    consider it an error.
    """
    if x.startswith('index file'):
        return True
    if x.startswith("WARNING"):
        return True
    return False


def _call_randomintersect(_self, other, iterations, intersect_kwargs,
                          shuffle_kwargs, report_iterations, debug,
                          _orig_processes):
    """
    Helper function that list-ifies the output from randomintersection, s.t.
    it can be pickled across a multiprocess Pool.
    """
    return list(
        _self.randomintersection(
            other, iterations,
            intersect_kwargs=intersect_kwargs,
            shuffle_kwargs=shuffle_kwargs,
            report_iterations=report_iterations,
            debug=False, processes=None,
            _orig_processes=_orig_processes)
    )


def close_or_delete(*args):
    """
    Single function that can be used to get rid of a BedTool, whether it's a
    streaming or file-based version.
    """
    for x in args:
        if isinstance(x.fn, six.string_types):
            os.unlink(x.fn)
        elif hasattr(x.fn, 'close'):
            x.fn.close()
        if hasattr(x.fn, 'throw'):
            x.fn.throw(StopIteration)


def n_open_fds():
    pid = os.getpid()
    procs = subprocess.check_output(
        ['lsof', '-w', '-Ff', '-p', str(pid)])
    nprocs = 0
    for i in procs.splitlines():
        if i[1:].isdigit() and i[0] == 'f':
            nprocs += 1
    return nprocs


import re
coord_re = re.compile(
    r"""
    (?P<chrom>.+):
    (?P<start>\d+)-
    (?P<stop>\d+)
    (?:\[(?P<strand>.)\])?""", re.VERBOSE)


def string_to_interval(s):
    """
    Convert string of the form "chrom:start-stop" or "chrom:start-stop[strand]"
    to an interval.

    Assumes zero-based coords.

    If it's already an interval, then return it as-is.
    """
    if isinstance(s, six.string_types):
        m = coord_re.search(s)
        if m.group('strand'):
            return create_interval_from_list([
                m.group('chrom'),
                m.group('start'),
                m.group('stop'),
                '.',
                '0',
                m.group('strand')])
        else:
            return create_interval_from_list([
                m.group('chrom'),
                m.group('start'),
                m.group('stop'),
            ])
    return s


def internet_on(timeout=1):
    try:
        response = urllib.request.urlopen('http://google.com', timeout=timeout)
        return True
    except urllib.error.URLError as err:
        pass
    return False


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
            print((subprocess.list2cmdline(cmds)))

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
    if isinstance(chrom_sizes, str):
        chrom_sizes = chromsizes(chrom_sizes)
    fout = open(fn, 'w')
    for chrom, bounds in sorted(chrom_sizes.items()):
        line = chrom + '\t' + str(bounds[1]) + '\n'
        fout.write(line)
    fout.close()
    return fn


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

atexit.register(cleanup)
