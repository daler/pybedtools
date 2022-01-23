import sys
import os
import gzip
import tempfile
import subprocess
import random
import string
import glob
import struct
import atexit
import six
import pysam
import re
from six.moves import urllib

try:  # Use genomepy to determine chrom sizes if it is installed
    import genomepy
except ImportError:
    pass
from . import cbedtools
from . import settings
from . import filenames
from . import genome_registry
from .logger import logger
from .cbedtools import create_interval_from_list

BUFSIZE = -1

_tags = {}


def set_bedtools_path(path=""):
    """
    Explicitly set path to `BEDTools` installation dir.

    If BEDTools is not available on your system path, specify the path to the
    dir containing the BEDTools executables (intersectBed, subtractBed, etc)
    with this function.

    To reset and use the default system path, call this function with no
    arguments or use path="".
    """
    from . import paths

    paths._set_bedtools_path(path)


def get_bedtools_path():
    """
    Returns the currently-set path to bedtools
    """
    from . import paths

    return paths._get_bedtools_path()


def set_R_path(path=""):
    """
    Explicitly set path to `R` installation dir.

    If R is not available on the path, then it can be explicitly
    specified here.

    Use path="" to reset to default system path.
    """
    from . import paths

    paths._set_R_path(path)


def _check_for_bedtools(force_check=False, verbose=False, override=None):
    """
    Checks installation as well as version (based on whether or not "bedtools
    intersect" works, or just "intersectBed")
    """
    if settings._bedtools_installed and not force_check:
        return True

    if (len(settings.bedtools_version) == 0) or force_check:
        try:
            v = (
                subprocess.check_output(
                    [os.path.join(settings._bedtools_path, "bedtools"), "--version"]
                )
                .decode("utf-8")
                .rstrip()
            )

            if verbose:
                print("Found bedtools version '%s'" % v)

            settings._bedtools_installed = True

            # Override, used for testing
            if override is not None:
                v = override

            # To allow more complicated versions as found in Linux distributions, e.g.:
            #  bedtools v2.26.0
            #  bedtools debian/2.28.0+dfsg-2-dirty
            m = re.search("^bedtools [^0-9]*([0-9][0-9.]*)", v)
            if not m:
                raise ValueError('Cannot identify version number from "{}"'.format(v))
            vv = m.group(1)

            settings.bedtools_version = [int(i) for i in vv.split(".")]

            settings._v_2_27_plus = (
                settings.bedtools_version[0] >= 2 and settings.bedtools_version[1] >= 27
            )

            settings._v_2_15_plus = (
                settings.bedtools_version[0] >= 2 and settings.bedtools_version[1] >= 15
            )

        except subprocess.CalledProcessError:
            if settings._bedtools_path:
                add_msg = "(tried path '%s')" % settings._bedtools_path
            else:
                add_msg = ""
            raise OSError(
                "Please make sure you have installed BEDTools"
                "(https://github.com/arq5x/bedtools) and that "
                "it's on the path. %s" % add_msg
            )


def _check_for_R():
    try:
        p = subprocess.Popen(
            [os.path.join(settings._R_path, "R"), "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        settings._R_installed = True
    except OSError:
        if settings._R_path:
            add_msg = "(tried path '%s')" % settings._R_path
        else:
            add_msg = ""
        raise ValueError("Please install R and ensure it is on your path %s" % add_msg)


class Error(Exception):
    """Base class for this module's exceptions"""

    pass


class pybedtoolsError(Error):
    pass


class BEDToolsError(Error):
    def __init__(self, cmd, msg):
        self.cmd = str(cmd)
        self.msg = str(msg)

    def __str__(self):
        m = (
            "\nCommand was:\n\n\t"
            + self.cmd
            + "\n"
            + "\nError message was:\n"
            + self.msg
        )
        return m


def isGZIP(fn):
    with open(fn, "rb") as f:
        start = f.read(3)
        if start == b"\x1f\x8b\x08":
            return True
    return False


def isBGZIP(fn):
    """
    Reads a filename to see if it's a BGZIPed file or not.
    """
    with open(fn, "rb") as fh:
        header_str = fh.read(15)

    if len(header_str) < 15:
        return False

    header = struct.unpack_from("BBBBiBBHBBB", header_str)

    id1, id2, cm, flg, mtime, xfl, os_, xlen, si1, si2, slen = header
    if (
        (id1 == 31)
        and (id2 == 139)
        and (cm == 8)
        and (flg == 4)
        and (si1 == 66)
        and (si2 == 67)
        and (slen == 2)
    ):
        return True
    return False


def isBAM(fn):
    """
    Returns True if the file is both BGZIPed and the compressed contents have
    start with the magic number `BAM\\x01`, or if the file is CRAM format (see
    isCRAM()).
    """
    # Note: previously we were catching ValueError when trying to open
    # a non-BAM with pysam.Samfile. That started segfaulting, so now do it the
    # right way with magic number.
    with gzip.open(fn, "rb") as in_:
        if isBGZIP(fn) and (in_.read(4).decode() == "BAM\x01"):
            return True
    if isCRAM(fn):
        return True


def isCRAM(fn):
    """
    Returns True if the file starts with the bytes for the characters "CRAM".
    """
    with open(fn, "rb") as in_:
        if in_.read(4).decode(errors="ignore") == "CRAM":
            return True


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
        errstr = "The tempdir you specified, %s, does not exist" % tempdir
        if six.PY2:
            raise ValueError(errstr)
        raise FileNotFoundError(errstr)
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
            print("removing", fn)
        if os.path.exists(fn):
            os.unlink(fn)
    if remove_all:
        fns = glob.glob(os.path.join(get_tempdir(), "pybedtools.*.tmp"))
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
        raise BEDToolsError(prog_name, prog_name + "not a recognized BEDTools program")
    return [os.path.join(settings._bedtools_path, "bedtools"), prog_name]


def call_bedtools(
    cmds,
    tmpfn=None,
    stdin=None,
    check_stderr=None,
    decode_output=True,
    encode_input=True,
):
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

    *decode_output* should be set to False when you are iterating over a BAM
    file, where the data represent binary rather than text data.
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
                "helpers.call_bedtools(): input is stream, output is " "stream"
            )
            logger.debug("helpers.call_bedtools(): cmds=%s", " ".join(cmds))
            p = subprocess.Popen(
                cmds,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                stdin=subprocess.PIPE,
                bufsize=BUFSIZE,
            )
            if encode_input:
                for line in stdin:
                    p.stdin.write(line.encode())
            else:
                for line in stdin:
                    p.stdin.write(line)

            # This is important to prevent deadlocks
            p.stdin.close()

            if decode_output:
                output = (i.decode("UTF-8") for i in p.stdout)
            else:
                output = (i for i in p.stdout)

            stderr = None

        # coming from an iterator, writing to file
        if input_is_stream and not output_is_stream:
            logger.debug("helpers.call_bedtools(): input is stream, output is file")
            logger.debug("helpers.call_bedtools(): cmds=%s", " ".join(cmds))
            outfile = open(tmpfn, "wb")
            p = subprocess.Popen(
                cmds,
                stdout=outfile,
                stderr=subprocess.PIPE,
                stdin=subprocess.PIPE,
                bufsize=BUFSIZE,
            )
            if hasattr(stdin, "read"):
                stdout, stderr = p.communicate(stdin.read())
            else:
                for item in stdin:
                    p.stdin.write(item.encode())
                stdout, stderr = p.communicate()
            output = tmpfn
            outfile.close()

        # coming from a file, sending as iterator
        if not input_is_stream and output_is_stream:
            logger.debug(
                "helpers.call_bedtools(): input is filename, " "output is stream"
            )
            logger.debug("helpers.call_bedtools(): cmds=%s", " ".join(cmds))
            p = subprocess.Popen(
                cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=BUFSIZE
            )
            if decode_output:
                output = (i.decode("UTF-8") for i in p.stdout)
            else:
                output = (i for i in p.stdout)
            stderr = None

        # file-to-file
        if not input_is_stream and not output_is_stream:
            logger.debug(
                "helpers.call_bedtools(): input is filename, output "
                "is filename (%s)",
                tmpfn,
            )
            cmds = list(map(str, cmds))
            logger.debug("helpers.call_bedtools(): cmds=%s", " ".join(cmds))
            outfile = open(tmpfn, "wb")
            p = subprocess.Popen(
                cmds, stdout=outfile, stderr=subprocess.PIPE, bufsize=BUFSIZE
            )
            stdout, stderr = p.communicate()
            output = tmpfn
            outfile.close()

        # Check if it's OK using a provided function to check stderr. If it's
        # OK, dump it to sys.stderr so it's printed, and reset it to None so we
        # don't raise an exception
        if check_stderr is not None:
            if isinstance(stderr, bytes):
                stderr = stderr.decode("UTF_8")
            if check_stderr(stderr):
                sys.stderr.write(stderr)
                stderr = None

        if stderr:
            # Fix for issue #147. In general, we consider warnings to not be
            # fatal, so just show 'em and continue on.
            #
            # bedtools source has several different ways of showing a warning,
            # but they seem to all have "WARNING" in the first 20 or so
            # characters
            if isinstance(stderr, bytes):
                stderr = stderr.decode("UTF_8")
            if len(stderr) > 20 and "WARNING" in stderr[:20].upper():
                sys.stderr.write(stderr)
            else:
                raise BEDToolsError(subprocess.list2cmdline(cmds), stderr)

    except (OSError, IOError) as err:
        print("%s: %s" % (type(err), os.strerror(err.errno)))
        print("The command was:\n\n\t%s\n" % subprocess.list2cmdline(cmds))

        problems = {
            2: (
                "* Did you spell the command correctly?",
                "* Do you have BEDTools installed and on the path?",
            ),
            13: (
                "* Do you have permission to write "
                'to the output file ("%s")?' % tmpfn,
            ),
            24: (
                "* Too many files open -- please submit "
                "a bug report so that this can be fixed",
            ),
            32: (
                "* Broken pipe -- if you passed a BedTool object "
                "that was created using a generator function, "
                "please try saving it to disk first using the "
                ".saveas() method before calling this bedtools "
                "command. See issue #49 for more.",
            ),
        }

        print("Things to check:")
        print("\n\t" + "\n\t".join(problems[err.errno]))
        raise OSError("See above for commands that gave the error")

    return output


def _check_sequence_stderr(x):
    """
    If stderr created by fastaFromBed starts with 'index file', then don't
    consider it an error.
    """
    if isinstance(x, bytes):
        x = x.decode("UTF-8")
    if x.startswith("index file"):
        return True
    if x.startswith("WARNING"):
        return True
    return False


def _call_randomintersect(
    _self,
    other,
    iterations,
    intersect_kwargs,
    shuffle_kwargs,
    report_iterations,
    debug,
    _orig_processes,
):
    """
    Helper function that list-ifies the output from randomintersection, s.t.
    it can be pickled across a multiprocess Pool.
    """
    return list(
        _self.randomintersection(
            other,
            iterations,
            intersect_kwargs=intersect_kwargs,
            shuffle_kwargs=shuffle_kwargs,
            report_iterations=report_iterations,
            debug=False,
            processes=None,
            _orig_processes=_orig_processes,
        )
    )


def close_or_delete(*args):
    """
    Single function that can be used to get rid of a BedTool, whether it's a
    streaming or file-based version.
    """
    for x in args:
        if isinstance(x.fn, six.string_types):
            os.unlink(x.fn)
        elif hasattr(x.fn, "close"):
            x.fn.close()
        if hasattr(x.fn, "throw"):
            x.fn.throw(StopIteration)


def n_open_fds():
    pid = os.getpid()
    procs = subprocess.check_output(["lsof", "-w", "-Ff", "-p", str(pid)])
    nprocs = 0
    for i in procs.splitlines():
        if i[1:].isdigit() and i[0] == "f":
            nprocs += 1
    return nprocs


import re

coord_re = re.compile(
    r"""
    (?P<chrom>.+):
    (?P<start>\d+)-
    (?P<stop>\d+)
    (?:\[(?P<strand>.)\])?""",
    re.VERBOSE,
)


def string_to_interval(s):
    """
    Convert string of the form "chrom:start-stop" or "chrom:start-stop[strand]"
    to an interval.

    Assumes zero-based coords.

    If it's already an interval, then return it as-is.
    """
    if isinstance(s, six.string_types):
        m = coord_re.search(s)
        if m.group("strand"):
            return create_interval_from_list(
                [
                    m.group("chrom"),
                    m.group("start"),
                    m.group("stop"),
                    ".",
                    "0",
                    m.group("strand"),
                ]
            )
        else:
            return create_interval_from_list(
                [m.group("chrom"), m.group("start"), m.group("stop")]
            )
    return s


class FisherOutput(object):
    def __init__(self, s, **kwargs):
        """
        fisher returns text results like::

            # Contingency Table
            #_________________________________________
            #           | not in -b    | in -b        |
            # not in -a | 3137160615   | 503          |
            #     in -a | 100          | 46           |
            #_________________________________________
            # p-values for fisher's exact test
            left	right	two-tail	ratio
            1.00000	0.00000	0.00000	2868973.922

        """
        if isinstance(s, str):
            s = open(s).read()
        if hasattr(s, "next"):
            s = "".join(i for i in s)
        table = {
            "not in -a": {"not in -b": None, "in -b": None},
            "in -a": {"not in -b": None, "in -b": None},
        }
        self.text = s
        lines = s.splitlines()
        for i in lines:
            if "not in -a" in i:
                _, in_b, not_in_b, _ = i.strip().split("|")
                table["not in -a"]["not in -b"] = int(not_in_b)
                table["not in -a"]["in -b"] = int(in_b)

            if "    in -a" in i:
                _, in_b, not_in_b, _ = i.strip().split("|")
                table["in -a"]["not in -b"] = int(not_in_b)
                table["in -a"]["in -b"] = int(in_b)
        self.table = table
        left, right, two_tail, ratio = lines[-1].split()
        self.left_tail = float(left)
        self.right_tail = float(right)
        self.two_tail = float(two_tail)
        self.ratio = float(ratio)

    def __str__(self):
        return self.text

    def __repr__(self):
        return "<%s at %s>\n%s" % (self.__class__.__name__, id(self), self.text)


class SplitOutput(object):
    def __init__(self, output, **kwargs):
        """
        Handles output from bedtools split, which sends a report of files to
        stdout. This class parses that list into something more convenient to
        use within pybedtools.

        Most useful is probably the .bedtools attribute, which is a list of
        BedTool objects.
        """
        from .bedtool import BedTool

        if isinstance(output, str):
            output = open(output).read()
        if hasattr(output, "next"):
            output = "".join(i for i in output)

        #: store a copy of the output
        self.text = output

        #: BedTool objects created from output
        self.bedtools = []

        #: Filenames that were created from the split
        self.files = []

        #: number of bases in each file
        self.nbases = []

        #: number of features in each file
        self.counts = []

        for line in output.splitlines():
            toks = line.split()
            self.files.append(toks[0])
            self.nbases.append(int(toks[1]))
            self.counts.append(int(toks[2]))
            self.bedtools.append(BedTool(toks[0]))


def internet_on(timeout=1):
    try:
        response = urllib.request.urlopen("http://genome.ucsc.edu", timeout=timeout)
        return True
    except urllib.error.URLError as err:
        pass
    return False


def get_chromsizes_from_ucsc(
    genome,
    saveas=None,
    mysql="mysql",
    fetchchromsizes="fetchChromSizes",
    timeout=None,
    host_url="genome-mysql.cse.ucsc.edu",
):
    """
    Download chrom size info for *genome* from UCSC and returns the dictionary.

    Parameters
    ----------

    genome : str
        Name of the genome assembly (e.g., "hg38")

    saveas : str
        Filename to save output to. Dictionary will still be returned.

    mysql, fetchchromsizes : str
        Paths to MySQL and fetchChromSizes.

    timeout : float
        How long to wait for a response; mostly used for testing.

    host_url : str
        URL of UCSC mirror MySQL server.
    """
    if not internet_on(timeout=timeout):
        raise ValueError(
            "It appears you don't have an internet connection "
            "-- unable to get chromsizes from UCSC"
        )
    cmds = [
        mysql,
        "--user=genome",
        "--host=" + host_url,
        "-A",
        "-e",
        "select chrom, size from %s.chromInfo" % genome,
    ]
    failures = []
    d = {}
    try:
        p = subprocess.Popen(
            cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=BUFSIZE
        )
        stdout, stderr = p.communicate()
        if stderr:
            print(stderr)
            print("Commands were:\n")
            print((subprocess.list2cmdline(cmds)))

        lines = stdout.splitlines()[1:]
        for line in lines:
            if isinstance(line, bytes):
                line = line.decode("UTF-8")
            chrom, size = line.split()
            d[chrom] = (0, int(size))

        if saveas is not None:
            chromsizes_to_file(d, saveas)

    except OSError as err:
        if err.errno == 2:
            failures.append("Can't find mysql at path {0}".format(mysql))
        else:
            raise
    try:
        cmds = [fetchchromsizes, genome]
        p = subprocess.Popen(
            cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=BUFSIZE
        )
        stdout, stderr = p.communicate()
        if stderr:
            if "INFO: trying WGET" not in str(stderr):
                print(stderr)
                print("Commands were:\n")
                print((subprocess.list2cmdline(cmds)))

        lines = stdout.splitlines()
        for line in lines:
            if isinstance(line, bytes):
                line = line.decode("UTF-8")
            chrom, size = line.split()
            d[chrom] = (0, int(size))

        if saveas is not None:
            chromsizes_to_file(d, saveas)

    except OSError as err:
        if err.errno == 2:
            failures.append("Can't find path to fetchChromsizes")

    if not d:
        raise OSError(failures)
    return d


def chromsizes_to_file(chrom_sizes, fn=None):
    """
    Converts a *chromsizes* dictionary to a file.  If *fn* is None, then a
    tempfile is created (which can be deleted with pybedtools.cleanup()).

    Returns the filename.
    """
    if fn is None:
        tmpfn = tempfile.NamedTemporaryFile(
            prefix="pybedtools.", suffix=".tmp", delete=False
        )
        tmpfn = tmpfn.name
        filenames.TEMPFILES.append(tmpfn)
        fn = tmpfn
    if isinstance(chrom_sizes, str):
        chrom_sizes = chromsizes(chrom_sizes)
    fout = open(fn, "wt")
    for chrom, bounds in sorted(chrom_sizes.items()):
        line = chrom + "\t" + str(bounds[1]) + "\n"
        fout.write(line)
    fout.close()
    return fn


def get_chromsizes_from_genomepy(
    genome, saveas=None,
):
    """
    Get chrom size info for *genome* from genomepy, if genomepy is installed.

    Parameters
    ----------

    genome : str
        Name of the genome assembly (e.g., "hg38")

    saveas : str
        Filename to save output to. Dictionary will still be returned.
    """
    if "genomepy" not in sys.modules:
        return None

    d = {}
    try:
        g = genomepy.Genome(genome)
        # Fail silently if the sizes file cannot be accessed
        if not hasattr(g, "sizes_file"):
            return None
        for line in open(g.sizes_file):
            chrom, size = line.split()
            d[chrom] = (0, int(size))

        if saveas is not None:
            chromsizes_to_file(d, saveas)
    except FileNotFoundError:
        return None

    return d


def chromsizes(genome):
    """
    Looks for a *genome* already included in the genome registry; if not found
    it first tries to look it up via genomepy. If genomepy is not installed, or
    if this lookup fails then it looks it up on UCSC.  Returns the dictionary of 
    chromsize tuples where each tuple has (start,stop).

    Chromsizes are described as (start, stop) tuples to allow randomization
    within specified regions; e. g., you can make a chromsizes dictionary that
    represents the extent of a tiling array.

    Example usage:

        >>> dm3_chromsizes = chromsizes('dm3')
        >>> for i in sorted(dm3_chromsizes.items()):
        ...     print(i)
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
        chromsizes = get_chromsizes_from_genomepy(genome)
        if chromsizes is None:
            return get_chromsizes_from_ucsc(genome)
        else:
            return chromsizes


def get_includes():
    """
    Returns a list of include directories with BEDTools headers
    """
    dirname = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    return [dirname, os.path.join(dirname, "include")]


atexit.register(cleanup)
