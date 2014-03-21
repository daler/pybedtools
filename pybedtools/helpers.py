import sys
import os
import tempfile
import subprocess
import random
import string
import glob
import struct
import atexit
import pybedtools

import settings

BUFSIZE = 1

_tags = {}


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
    for key, item in _tags.iteritems():
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


class History(list):
    def __init__(self):
        """
        Represents one or many HistorySteps.  Mostly used for nicely formatting
        a series of HistorySteps.
        """
        list.__init__(self)


class HistoryStep(object):
    def __init__(self, method, args, kwargs, bedtool_instance,
                 parent_tag, result_tag):
        """
        Class to represent one step in the history.

        Mostly used for its __repr__ method, to try and exactly replicate code
        that can be pasted to re-do history steps
        """
        try:
            self.method = method._name
        except AttributeError:
            self.method = method.func_name
        self.args = args
        self.kwargs = kwargs
        self.fn = bedtool_instance.fn
        tag = ''.join(random.choice(string.lowercase) for _ in xrange(8))
        self.parent_tag = parent_tag
        self.result_tag = result_tag

    def _clean_arg(self, arg):
        """
        Wrap strings in quotes and convert bedtool instances to filenames.
        """
        if isinstance(arg, pybedtools.BedTool):
            arg = arg.fn
        if isinstance(arg, basestring):
            arg = '"%s"' % arg
        return arg

    def __repr__(self):
        # Still not sure whether to use pybedtools.bedtool() or bedtool()
        s = ''
        s += '<HistoryStep> '
        if os.path.exists(self.fn):
            s += 'BedTool("%(fn)s").%(method)s(%%s%%s)' % self.__dict__
        else:
            s += 'BedTool("MISSING FILE: %(fn)s")' % self.__dict__
            s += '.%(method)s(%%s%%s)' % self.__dict__

        # Format args and kwargs
        args_string = ','.join(map(self._clean_arg, self.args))
        kwargs_string = ','.join(
            ['%s=%s' % (i[0], self._clean_arg(i[1]))
             for i in self.kwargs.items()])
        # stick a comma on the end if there's something here
        if len(args_string) > 0:
            args_string += ', '

        s = s % (args_string, kwargs_string)
        s += ', parent tag: %s' % self.parent_tag
        s += ', result tag: %s' % self.result_tag
        return s


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
    for fn in pybedtools.BedTool.TEMPFILES:
        if verbose:
            print 'removing', fn
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
            pybedtools.logger.debug(
                'helpers.call_bedtools(): input is stream, output is '
                'stream')
            pybedtools.logger.debug(
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
            pybedtools.logger.debug(
                'helpers.call_bedtools(): input is stream, output is file')
            pybedtools.logger.debug(
                'helpers.call_bedtools(): cmds=%s', ' '.join(cmds))
            outfile = open(tmpfn, 'w')
            p = subprocess.Popen(cmds,
                                 stdout=outfile,
                                 stderr=subprocess.PIPE,
                                 stdin=subprocess.PIPE,
                                 bufsize=BUFSIZE)
            if isinstance(stdin, file):
                stdout, stderr = p.communicate(stdin.read())
            else:
                for item in stdin:
                    p.stdin.write(item)
                stdout, stderr = p.communicate()
            output = tmpfn
            outfile.close()

        # coming from a file, sending as iterator
        if not input_is_stream and output_is_stream:
            pybedtools.logger.debug(
                'helpers.call_bedtools(): input is filename, '
                'output is stream')
            pybedtools.logger.debug(
                'helpers.call_bedtools(): cmds=%s', ' '.join(cmds))
            p = subprocess.Popen(cmds,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 bufsize=BUFSIZE)
            output = p.stdout
            stderr = None

        # file-to-file
        if not input_is_stream and not output_is_stream:
            pybedtools.logger.debug(
                'helpers.call_bedtools(): input is filename, output '
                'is filename (%s)', tmpfn)
            pybedtools.logger.debug(
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
        print '%s: %s' % (type(err), os.strerror(err.errno))
        print 'The command was:\n\n\t%s\n' % subprocess.list2cmdline(cmds)

        problems = {
            2: ('* Did you spell the command correctly?',
                '* Do you have BEDTools installed and on the path?'),
            13: ('* Do you have permission to write '
                 'to the output file ("%s")?' % tmpfn,),
            24: ('* Too many files open -- please submit '
                 'a bug report so that this can be fixed',)
        }

        print 'Things to check:'
        print '\n\t' + '\n\t'.join(problems[err.errno])
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
        if isinstance(x.fn, basestring):
            os.unlink(x.fn)
        elif hasattr(x.fn, 'close'):
            x.fn.close()
        if hasattr(x.fn, 'throw'):
            x.fn.throw(StopIteration)

def _jaccard_output_to_dict(s, **kwargs):
    """
    jaccard method doesn't return an interval file, rather, it returns a short
    summary of results.  Here, we simply parse it into a dict for convenience.
    """
    if isinstance(s, basestring):
        s = open(s).read()
    if hasattr(s, 'next'):
        s = ''.join(i for i in s)
    header, data = s.splitlines()
    header = header.split()
    data = data.split()
    data[0] = int(data[0])
    data[1] = int(data[1])
    data[2] = float(data[2])
    data[3] = int(data[3])
    return dict(zip(header, data))


def _reldist_output_handler(s, **kwargs):
    """
    reldist, if called with -detail, returns a valid BED file with the relative
    distance as the last field.  In that case, return the BedTool immediately.
    If not -detail, then the results are a table, in which case here we parse
    into a dict for convenience.
    """
    if 'detail' in kwargs:
        return pybedtools.BedTool(s)
    if isinstance(s, basestring):
        iterable = open(s)
    if hasattr(s, 'next'):
        iterable = s
    header = iterable.next().split()
    results = {}
    for h in header:
        results[h] = []
    for i in iterable:
        reldist, count, total, fraction = i.split()
        data = [
            float(reldist),
            int(count),
            int(total),
            float(fraction)
        ]
        for h, d in zip(header, data):
            results[h].append(d)
    return results


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
    if isinstance(s, basestring):
        m = coord_re.search(s)
        if m.group('strand'):
            return pybedtools.create_interval_from_list([
                m.group('chrom'),
                m.group('start'),
                m.group('stop'),
                '.',
                '0',
                m.group('strand')])
        else:
            return pybedtools.create_interval_from_list([
                m.group('chrom'),
                m.group('start'),
                m.group('stop'),
            ])
    return s

atexit.register(cleanup)
