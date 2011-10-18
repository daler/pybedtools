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

# Check calls against these names to only allow calls to known BEDTools
# programs (basic security)
_prog_names = ['annotateBed', 'bedToBam', 'complementBed', 'flankBed',
'linksBed', 'overlap', 'shuffleBed', 'subtractBed', 'bamToBed', 'bedToIgv',
'coverageBed', 'genomeCoverageBed', 'maskFastaFromBed', 'pairToBed', 'slopBed',
'unionBedGraphs', 'bed12ToBed6', 'closestBed', 'fastaFromBed', 'intersectBed',
'mergeBed', 'pairToPair', 'sortBed', 'windowBed', 'groupBy', 'tagBam',
'nucBed', 'multiBamCov']

_tags = {}


class Error(Exception):
    """Base class for this module's exceptions"""
    pass


class BEDToolsError(Error):
    pass


def isBAM(fn):
    """
    Reads a filename to see if it's a BAM file or not.
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
        kwargs_string = ','.join(['%s=%s' % \
                        (i[0], self._clean_arg(i[1])) \
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
    Sets the directory for temp files.  Useful for clusters that use a /scratch
    partition rather than a /tmp dir.  Convenience function to simply set
    tempfile.tempdir.
    """
    if not os.path.exists(tempdir):
        errstr = 'The tempdir you specified, %s, does not exist' % tempdir
        raise ValueError(errstr)
    tempfile.tempdir = tempdir


def get_tempdir():
    """
    Gets the current tempdir for the module.
    """
    return tempfile.tempdir


def cleanup(verbose=False, remove_all=False):
    """
    Deletes all temporary files in the *BedTool.TEMPFILES* class
    variable.

    If *verbose*, reports what it's doing

    If *remove_all*, then ALL files matching "pybedtools.*.tmp" in the temp dir
    will be deleted.
    """
    if pybedtools.KEEP_TEMPFILES:
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
    instream = stdin is not None
    outstream = tmpfn is None

    if cmds[0] not in _prog_names:
        raise BEDToolsError('"%s" not a recognized BEDTools program' % cmds[0])

    # use specifed path, "" by default. Special-case groupBy.
    if cmds[0] == 'groupBy':
        cmds[0] = os.path.join(pybedtools._filo_path, cmds[0])
    else:
        cmds[0] = os.path.join(pybedtools._path, cmds[0])

    try:
        # coming from an iterator, sending as iterator
        if instream and outstream:

            p = subprocess.Popen(cmds,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 stdin=subprocess.PIPE,
                                 bufsize=1)
            for line in stdin:
                p.stdin.write(line + "\n")

            # This is important to prevent deadlocks
            p.stdin.close()

            output = p.stdout
            stderr = None

        # coming from an iterator, writing to file
        if instream and not outstream:
            outfile = open(tmpfn, 'w')
            p = subprocess.Popen(cmds,
                                 stdout=outfile,
                                 stderr=subprocess.PIPE,
                                 stdin=subprocess.PIPE,
                                 bufsize=1)
            if isinstance(stdin, file):
                stdout, stderr = p.communicate(stdin.read())
            else:
                for item in stdin:
                    p.stdin.write(item + "\n")
                stdout, stderr = p.communicate()
            output = tmpfn
            outfile.close()

        # coming from a file, sending as iterator
        if not instream and outstream:
            p = subprocess.Popen(cmds,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 bufsize=1)
            output = p.stdout
            stderr = None

        # file-to-file
        if not instream and not outstream:
            outfile = open(tmpfn, 'w')
            p = subprocess.Popen(cmds,
                                 stdout=outfile,
                                 stderr=subprocess.PIPE,
                                 bufsize=1)
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
            sys.stderr.write('\nCommand was:\n\n\t%s\n' % \
                             subprocess.list2cmdline(cmds))
            sys.stderr.write('\nError message was:\n')
            sys.stderr.write(stderr)
            raise BEDToolsError('Error message from BEDTools written to '
                                'stderr, above', stderr)

    except (OSError, IOError) as err:
        print '%s: %s' % (type(err), os.strerror(err.errno))
        print 'The command was:\n\n\t%s\n' % subprocess.list2cmdline(cmds)

        problems = {2: ('* Did you spell the command correctly?',
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
    If BEDTools is not available on your system path, specify the path to the
    dir containing the BEDTools executables (intersectBed, subtractBed, etc)
    with this function.

    To reset and use the default system path, call this function with no
    arguments or use path="".
    """
    pybedtools._path = path


def set_samtools_path(path=""):
    """
    If samtools is not available on the path, then it can be explicitly
    specified here.

    Use path="" to reset to default system path.
    """
    pybedtools._samtools_path = path


def set_filo_path(path=""):
    """
    If filo is not available on the path, then it can be explicitly
    specified here.

    Use path="" to reset to default system path.
    """
    pybedtools._filo_path = path


def _check_sequence_stderr(x):
    """
    If stderr created by fastaFromBed starst with 'index file', then don't
    consider it an error.
    """
    if x.startswith('index file'):
        return True
    return False


def _call_randomintersect(_self, other, iterations, intersect_kwargs,
        shuffle_kwargs, report_iterations, debug, _orig_processes):
    """
    Helper function that list-ifies the output from randomintersection, s.t.
    it can be pickled across a multiprocess Pool.
    """
    return list(_self.randomintersection(other, iterations,
                                        intersect_kwargs=intersect_kwargs,
                                        shuffle_kwargs=shuffle_kwargs,
                                        report_iterations=report_iterations,
                                        debug=False, processes=None,
                                        _orig_processes=_orig_processes))

atexit.register(cleanup)
