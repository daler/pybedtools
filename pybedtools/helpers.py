import sys
import os
import tempfile
import subprocess
import random
import string
import glob
import pybedtools

# Check calls against these names to only allow calls to known BEDTools
# programs (basic security)
_prog_names = ['annotateBed', 'bedToBam', 'complementBed', 'flankBed',
'linksBed', 'overlap', 'shuffleBed', 'subtractBed', 'bamToBed', 'bedToIgv',
'coverageBed', 'genomeCoverageBed','maskFastaFromBed', 'pairToBed', 'slopBed',
'unionBedGraphs', 'bed12ToBed6', 'closestBed', 'fastaFromBed', 'intersectBed',
'mergeBed', 'pairToPair', 'sortBed', 'windowBed', ]

_tags = {}

class Error(Exception):
    """Base class for this module's exceptions"""
    pass

class BEDToolsError(Error):
    pass

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
    return '%s not found' % tag


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
        self.method = method.func_name
        self.args = args
        self.kwargs = kwargs
        self.fn = bedtool_instance.fn
        tag = ''.join(random.choice(string.lowercase) for _ in xrange(8))
        self.parent_tag = parent_tag
        self.result_tag = result_tag

    def _clean_arg(self,arg):
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
            s += 'BedTool("MISSING FILE: %(fn)s").%(method)s(%%s%%s)' % self.__dict__

        # Format args and kwargs
        args_string = ','.join(map(self._clean_arg, self.args))
        kwargs_string = ','.join(['%s=%s'% (i[0], self._clean_arg(i[1])) for i in self.kwargs.items()])

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
        raise ValueError, 'The tempdir you specified, %s, does not exist' % tempdir
    tempfile.tempdir = tempdir

def get_tempdir():
    """
    Gets the current tempdir for the module.
    """
    return tempfile.tempdir

def cleanup(verbose=True, remove_all=False):
    """
    Deletes all temporary files in the *BedTool.TEMPFILES* class
    variable.

    If *verbose*, reports what it's doing

    If *remove_all*, then ALL files matching "pybedtools.*.tmp" in the temp dir
    will be deleted.
    """
    for fn in pybedtools.BedTool.TEMPFILES:
        if verbose:
            print 'removing', fn
        if os.path.exists(fn):
            os.unlink(fn)
    if remove_all:
        fns = glob.glob(os.path.join(get_tempdir(), 'pybedtools.*.tmp'))
        for fn in fns:
            os.unlink(fn)


def _file_or_bedtool():
    '''
    Decorator that adds a line to the docstring indicating
    that a bedtool object is returned.
    '''
    extra_help = """
    .. note::

        This method accepts either a bedtool or a file name as the first
        unnamed argument

    """

    def decorator(func):
        """
        Adds the help to the function's __doc__
        """
        if func.__doc__ is None:
            func.__doc__ = ''
        orig = func.__doc__
        func.__doc__ += extra_help
        return func

    return decorator

def _returns_bedtool():
    '''
    Decorator that adds a line to the docstring indicating
    that a bedtool object is returned.
    '''
    extra_help = """
    .. note::

        This method returns a new bedtool instance
    """

    def decorator(func):
        """
        Adds the help to the function's __doc__
        """
        if func.__doc__ is None:
            func.__doc__ = ''
        orig = func.__doc__
        func.__doc__ += extra_help
        return func

    return decorator

def _implicit(option):
    '''
    Decorator that adds a line to the docstring indicating
    that a particular option is implied to be the default
    '''
    extra_help = """
    .. note::

        For convenience, the file this bedtool object points to is passed as "%s"
    """ % option

    def decorator(func):
        """
        Adds the help to the function's __doc__
        """
        if func.__doc__ is None:
            func.__doc__ = ''
        orig = func.__doc__
        func.__doc__ += extra_help
        return func

    return decorator

def _help(command):
    '''Decorator that adds help from each of the BEDtools programs to the
    docstring of the method that calls the program'''

    p = subprocess.Popen([command,'-h'], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    help_str = p.communicate()[1]
    help_str = help_str.replace('_','**')

    # insert tabs into the help
    help_str = help_str.split('\n')
    help_str = ['\t'+i for i in help_str]
    help_str = '\n'.join(help_str)

    def decorator(func):
        """
        Adds the help to the function's __doc__
        """
        if func.__doc__ is None:
            func.__doc__ = ''
        orig = func.__doc__
        func.__doc__ = '*pybedtools help:*\n'
        func.__doc__ += orig
        func.__doc__ += '\n\n*Original BEDtools program help:*\n'
        func.__doc__ += help_str
        return func

    return decorator

def call_bedtools(cmds, tmpfn=None, check_stderr=None):
    """
    Use subprocess.Popen to call BEDTools and catch any errors.

    Output goes to *tmpfn*, or, if None, output stays in subprocess.PIPE and
    can be iterated over.

    Prints some useful help upon getting common errors.

    *check_stderr* is a function that takes the stderr string as input and
    returns True if it's OK (that is, it's not really an error).  This is
    needed, e.g., for calling fastaFromBed which will report that it has to
    make a .fai for a fasta file.
    """
    if cmds[0] not in _prog_names:
        raise BEDToolsError('"%s" not a recognized BEDTools program' % cmds[0])
    try:
        if tmpfn is None:
            output = subprocess.PIPE
        else:
            output = open(tmpfn, 'w')
        p = subprocess.Popen(cmds, stdout=output, stderr=subprocess.PIPE)
        stdout,stderr = p.communicate()

        # Check if it's OK; if so dump it to sys.stderr and reset it to None so
        # we don't raise an exception
        if check_stderr is not None:
            if check_stderr(stderr):
                sys.stderr.write(stderr)
                stderr = None

        if stderr:
            print 'Command was:\n\n\t%s\n' % subprocess.list2cmdline(cmds)
            print 'Error message was:\n'
            #print '\n'.join([i for i in stderr.splitlines() if i.startswith('***')])
            print stderr
            raise BEDToolsError('See above for commands and error message', stderr)

    except (OSError, IOError) as err:
        print '%s: %s' % (type(err), os.strerror(err.errno))
        print 'The command was:\n\n\t%s\n' % subprocess.list2cmdline(cmds)

        problems = {2 :('* Did you spell the command correctly?', '* Do you have BEDTools installed and on the path?'),
                    13:('* Do you have permission to write to the output file ("%s")?' % tmpfn,),
                   }

        print 'Things to check:'
        print '\n\t'+'\n\t'.join(problems[err.errno])
        raise OSError('See above for commands that gave the error')

    if tmpfn is None:
        return output
    else:
        return tmpfn
