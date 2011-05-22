import tempfile
from math import floor, ceil
import os
import sys
import random
import string
from itertools import groupby

from pybedtools.helpers import _file_or_bedtool, _help, _implicit,\
    _returns_bedtool, get_tempdir, _tags,\
    History, HistoryStep, call_bedtools, _flatten_list, IntervalIterator

from cbedtools import IntervalFile
import pybedtools


class BedTool(object):
    TEMPFILES = []

    def __init__(self, fn, from_string=False):
        """
        Wrapper around Aaron Quinlan's ``BEDtools`` suite of programs
        (https://github.com/arq5x/bedtools); also contains many useful
        methods for more detailed work with BED files.

        *fn* is typically the name of a BED-like file, but can also be
        one of the following:

            * a string filename
            * another BedTool object
            * an iterable of Interval objects
            * an open file object
            * a "file contents" string (see below)

        If *from_string* is True, then you can pass a string that contains
        the contents of the BedTool you want to create.  This will treat all
        spaces as TABs and write to tempfile, treating whatever you pass as
        *fn* as the contents of the bed file.  This also strips empty lines.

        Typical usage is to point to an existing file::

            a = BedTool('a.bed')

        But you can also create one from scratch from a string::

            >>> s = '''
            ... chrX  1  100
            ... chrX 25  800
            ... '''
            >>> a = BedTool(s,from_string=True).saveas('a.bed')

        Or use examples that come with pybedtools::

             >>> example_files = pybedtools.list_example_files()
             >>> assert example_files[0] == 'a.bed'
             >>> a = pybedtools.example_bedtool('a.bed')

        """
        if not from_string:
            if isinstance(fn, BedTool):
                fn = fn.fn
            elif isinstance(fn, basestring):
                if not os.path.exists(fn):
                    raise ValueError('File "%s" does not exist' % fn)
            else:
                fn = fn
        else:
            bed_contents = fn
            fn = self._tmp()
            fout = open(fn, 'w')
            for line in bed_contents.splitlines():
                if len(line.strip()) == 0:
                    continue
                line = '\t'.join(line.split()) + '\n'
                fout.write(line)
            fout.close()

        tag = ''.join([random.choice(string.lowercase) for _ in xrange(8)])
        self._tag = tag
        _tags[tag] = self
        self.fn = fn
        self._hascounts = False

        self.history = History()

    def delete_temporary_history(self, ask=True, raw_input_func=None):
        """
        Use at your own risk!  This method will delete temp files. You will be
        prompted for deletion of files unless you specify *ask=False*.

        Deletes all temporary files created during the history of this BedTool
        up to but not including the file this current BedTool points to.

        Any filenames that are in the history and have the following pattern
        will be deleted::

            <TEMP_DIR>/pybedtools.*.tmp

        (where <TEMP_DIR> is the result from get_tempdir() and is by default
        "/tmp")

        Any files that don't have this format will be left alone.

        (*raw_input_func* is used for testing)
        """
        flattened_history = _flatten_list(self.history)
        to_delete = []
        tempdir = get_tempdir()
        for i in flattened_history:
            fn = i.fn
            if fn.startswith(os.path.join(os.path.abspath(tempdir),
                                          'pybedtools')):
                if fn.endswith('.tmp'):
                    to_delete.append(fn)

        if raw_input_func is None:
            raw_input_func = raw_input

        str_fns = '\n\t'.join(to_delete)
        if ask:
            answer = raw_input_func('Delete these files?\n\t%s\n(y/N) ' \
                                    % str_fns)

            if not answer.lower()[0] == 'y':
                print('OK, not deleting.')
                return
        for fn in to_delete:
            os.unlink(fn)
        return

    def _log_to_history(method):
        """
        Decorator to add a method and its kwargs to the history.

        Assumes that you only add this decorator to bedtool instances that
        return other bedtool instances
        """
        def decorated(self, *args, **kwargs):

            # this calls the actual method in the first place; *result* is
            # whatever you get back
            result = method(self, *args, **kwargs)

            # add appropriate tags
            parent_tag = self._tag
            result_tag = result._tag

            # log the sucka
            history_step = HistoryStep(method, args, kwargs, self, parent_tag,
                                       result_tag)

            # only add the current history to the new bedtool if there's
            # something to add
            if len(self.history) > 0:
                result.history.append(self.history)

            # but either way, add this history step to the result.
            result.history.append(history_step)

            return result

        decorated.__doc__ = method.__doc__
        return decorated

    def filter(self, func, *args, **kwargs):
        """
        Takes a function *func* that is called for each feature
        in the `BedTool` object and returns only those
        for which the function returns True.

        *args and **kwargs are passed directly to *func*.

        Returns a streaming BedTool; if you want the filename then use the
        .saveas() method.

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> subset = a.filter(lambda b: b.chrom == 'chr1' and b.start < 150)
        >>> len(a), len(subset)
        (4, 2)

        so it has extracted 2 records from the original 4.

        """
        return BedTool((f for f in self if func(f, *args, **kwargs)))

    def field_count(self, n=10):
        """
        Return the number of fields in the features this file contains.  Checks
        the first *n* features.
        """
        i = 0
        fields = set([])
        for feat in self:
            if i > n:
                break
            i += 1
            # TODO: make this more efficient.
            fields.update([len(feat.fields)])
        assert len(fields) == 1, fields
        return list(fields)[0]

    def each(self, func, *args, **kwargs):
        """
        Applies user-defined function *func* to each feature.  *func* must
        accept an Interval as its first argument; *args and **kwargs will be
        passed to *func*.

        *func* must return an Interval object.

        >>> def truncate_feature(feature, limit=0):
        ...     feature.score = str(len(feature))
        ...     if len(feature) > limit:
        ...         feature.stop = feature.start + limit
        ...         feature.name = feature.name + '.short'
        ...     return feature

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> b = a.each(truncate_feature, limit=100)
        >>> print b #doctest: +NORMALIZE_WHITESPACE
        chr1	1	100	feature1	99	+
        chr1	100	200	feature2	100	+
        chr1	150	250	feature3.short	350	-
        chr1	900	950	feature4	50	+
        <BLANKLINE>

        """
        return BedTool((func(f, *args, **kwargs) for f in self))

    @_help('bed12ToBed6')
    @_file_or_bedtool()
    @_implicit('-i')
    @_returns_bedtool()
    @_log_to_history
    def bed6(self, **kwargs):
        """
        convert a BED12 to a BED6 file
        """
        if not 'i' in kwargs:
            kwargs['i'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='bed12ToBed6', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    @property
    def file_type(self):
        """
        Return the type of the current file.  One of ('bed','vcf','gff').

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> a.file_type
        'bed'
        """

        v = IntervalFile(self.fn)
        return v.file_type

    def cut(self, indexes):
        """
        Similar to unix `cut` except indexes are 0-based, must be a list
        and the columns are returned in the order requested.

        In addition, indexes can contain keys of the GFF/GTF attributes,
        in which case the values are returned. e.g. 'gene_name' will return the
        corresponding name from a GTF, or 'start' will return the start
        attribute of a BED Interval.

        See .with_column() if you need to do more complex operations.
        """
        fh = open(self._tmp(), "w")
        for f in self:
            print >>fh, "\t".join(map(str, [f[attr] for attr in indexes]))
        fh.close()
        return BedTool(fh.name)

    @classmethod
    def _tmp(self):
        '''
        Makes a tempfile and registers it in the BedTool.TEMPFILES class
        variable.  Adds a "pybedtools." prefix and ".tmp" extension for easy
        deletion if you forget to call pybedtools.cleanup().
        '''
        tmpfn = tempfile.NamedTemporaryFile(prefix='pybedtools.',
                                            suffix='.tmp', delete=False)
        tmpfn = tmpfn.name
        BedTool.TEMPFILES.append(tmpfn)
        return tmpfn

    def __iter__(self):
        """
        Dispatches the right iterator depending on how this BedTool was
        created
        """
        # Plain ol' filename
        if isinstance(self.fn, basestring):
            return IntervalFile(self.fn)

        # Open file, like subprocess.PIPE.
        if isinstance(self.fn, file):
            return IntervalIterator(self.fn)

        # Otherwise assume fn is already an iterable
        else:
            return self.fn

    def __repr__(self):
        if isinstance(self.fn, file):
            return '<BedTool(stream)>'
        if isinstance(self.fn, basestring):
            if os.path.exists(self.fn):
                return '<BedTool(%s)>' % self.fn
            else:
                return '<BedTool(MISSING FILE: %s)>' % self.fn
        else:
            return repr(self.fn)

    def __str__(self):
        """
        Different methods to return the string, depending on how the BedTool
        was created.  If self.fn is anything but a basestring, the iterable
        will be consumed.
        """
        if isinstance(self.fn, basestring):
            f = open(self.fn)
            s = f.read()
            f.close()
            return s
        elif isinstance(self.fn, file):
            return self.fn.read()
        else:
            return '\n'.join(str(i) for i in iter(self)) + '\n'

    def __len__(self):
        return self.count()

    def __eq__(self, other):
        if isinstance(other, basestring):
            other_str = other
        elif isinstance(other, BedTool):
            if not isinstance(self.fn, basestring) or not \
                            isinstance(other.fn, basestring):
                raise NotImplementedError('Testing equality only supported for'
                                          ' BedTools that point to files')
            other_str = open(other.fn).read()
        if open(self.fn).read() == other_str:
            return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    @_file_or_bedtool()
    def __add__(self, other):
        return self.intersect(other, u=True)

    @_file_or_bedtool()
    def __sub__(self, other):
        return self.intersect(other, v=True)

    def head(self, n=10):
        """
        Prints the first *n* lines
        """
        for i, line in enumerate(open(self.fn)):
            if i == (n):
                break
            print line,

    @_returns_bedtool()
    def set_chromsizes(self, chromsizes):
        """
        Set the chromsizes for this genome. If *chromsizes* is a string, it
        will be considered a genome assembly name.  If that assembly name is
        not available in pybedtools.genome_registry, then it will be searched
        for on the UCSC Genome Browser.

        Example usage:

            >>> hg19 = pybedtools.chromsizes('hg19')
            >>> a = pybedtools.example_bedtool('a.bed')
            >>> a = a.set_chromsizes(hg19)
            >>> print a.chromsizes['chr1']
            (0, 249250621)

        """
        if isinstance(chromsizes, basestring):
            self.chromsizes = pybedtools.chromsizes(chromsizes)
        elif isinstance(chromsizes, dict):
            self.chromsizes = chromsizes
        else:
            raise ValueError("Need to specify chromsizes either as a string"
                             " (assembly name) or a dictionary")
        return self

    def handle_kwargs(self, prog, **kwargs):
        """
        Handle most cases of BEDTool program calls, but leave the specifics
        up to individual methods.

        *prog* is a BEDTools program name, e.g., 'intersectBed'.

        *kwargs* are passed directly from the calling method (like
        self.intersect).

        This method figures out, given how this BedTool was constructed, what
        to send to BEDTools programs -- for example, an open file to stdin with
        the `-` argument, or a filename with the `-a` argument.
        """
        # Dict of programs and which arguments *self.fn* can be used as
        implicit_instream1 = {'intersectBed': 'a',
                               'subtractBed': 'a',
                                'closestBed': 'a',
                                 'windowBed': 'a',
                                   'slopBed': 'i',
                                  'mergeBed': 'i',
                                   'sortBed': 'i',
                               'bed12ToBed6': 'i',
                                'shuffleBed': 'i',
                               'annotateBed': 'i',
                                  'flankBed': 'i',
                              'fastaFromBed': 'bed',
                          'maskFastaFromBed': 'bed',
                               'coverageBed': 'a'}

        # Which arguments *other.fn* can be used as
        implicit_instream2 = {'intersectBed': 'b',
                               'subtractBed': 'b',
                                'closestBed': 'b',
                                 'windowBed': 'b',
                               'coverageBed': 'b'}

        stdin = None

        # -----------------------------------------------------------------
        # Decide how to send instream1 to BEDTools.  If there's no implicit
        # instream1 arg, then do nothing.
        #
        try:
            # e.g., 'a' for intersectBed
            inarg1 = implicit_instream1[prog]

            # e.g., self.fn or 'a.bed' or an iterator...
            instream1 = kwargs[inarg1]

            # If it's a BedTool, then get underlying stream
            if isinstance(instream1, BedTool):
                instream1 = instream1.fn

            # Filename? No pipe, just provide the file
            if isinstance(instream1, basestring):
                kwargs[inarg1] = instream1
                stdin = None

            # Open file? Pipe it
            elif isinstance(instream1, file):
                kwargs[inarg1] = 'stdin'
                stdin = instream1

            # A generator or iterator: pipe it as a generator of lines
            else:
                kwargs[inarg1] = 'stdin'
                stdin = (str(i) for i in instream1)
        except KeyError:
            pass

        # -----------------------------------------------------------------
        # Decide how to send instream2 to BEDTools.
        try:
            # e.g., 'b' for intersectBed
            inarg2 = implicit_instream2[prog]

            # e.g., another BedTool
            instream2 = kwargs[implicit_instream2[prog]]

            # Get stream if BedTool
            if isinstance(instream2, BedTool):
                instream2 = instream2.fn

            # Filename
            if isinstance(instream2, basestring):
                kwargs[inarg2] = instream2

            # Otherwise we need to collapse it in order to send to BEDTools
            # programs
            else:
                collapsed_fn = self._tmp()
                fout = open(collapsed_fn, 'w')
                for i in instream2:
                    # TODO: does this need newlines?
                    fout.write(str(i))
                fout.close()
                kwargs[inarg2] = collapsed_fn
        except KeyError:
            pass

        # If stream not specified, then a tempfile will be created
        try:
            if kwargs.pop('stream'):
                tmp = None
            else:
                tmp = self._tmp()
        except KeyError:
            tmp = self._tmp()

        # Parse the kwargs into BEDTools-ready args
        cmds = [prog]
        for key, value in kwargs.items():
            if isinstance(value, bool):
                if value:
                    cmds.append('-' + key)
                else:
                    continue
            else:
                cmds.append('-' + key)
                cmds.append(str(value))
        return cmds, tmp, stdin

    @_returns_bedtool()
    @_log_to_history
    def remove_invalid(self):
        """
        Remove invalid features and return a new BedTool.

        >>> a = pybedtools.BedTool("chr1 10 100\\nchr1 10 1",
        ... from_string=True)
        >>> print a.remove_invalid() #doctest: +NORMALIZE_WHITESPACE
        chr1	10	100
        <BLANKLINE>

        """
        tmp = self._tmp()
        fout = open(tmp, 'w')
        i = iter(self)
        while True:
            try:
                fout.write(str(i.next()) + '\n')
            except pybedtools.MalformedBedLineError:
                continue
            except StopIteration:
                break
        fout.close()
        return BedTool(tmp)

    @_help('intersectBed')
    @_file_or_bedtool()
    @_implicit('-a')
    @_returns_bedtool()
    @_log_to_history
    def intersect(self, b=None, **kwargs):
        """
        Intersect with another BED file. If you want to use BAM as input, you
        need to specify *abam='filename.bam'*.  Returns a new BedTool object.

        Example usage:

            Create new BedTool object

            >>> a = pybedtools.example_bedtool('a.bed')

            Get overlaps with `b.bed`:

            >>> b = pybedtools.example_bedtool('b.bed')
            >>> overlaps = a.intersect(b)

            Use `v=True` to get the inverse -- those unique to "a.bed":

            >>> unique_to_a = a.intersect(b, v=True)
        """
        kwargs['b'] = b

        if ('abam' not in kwargs) and ('a' not in kwargs):
            kwargs['a'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='intersectBed', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    @_help('fastaFromBed')
    @_implicit('-bed')
    @_returns_bedtool()
    def sequence(self, fi, **kwargs):
        '''
        Wraps ``fastaFromBed``.  *fi* is passed in by the user; *bed* is
        automatically passed in as the bedfile of this object; *fo* by default
        is a temp file.  Use save_seqs() to save as a file.

        The end result is that this BedTool will have an attribute, self.seqfn,
        that points to the new fasta file.

        Example usage:

        >>> a = pybedtools.BedTool("""
        ... chr1 1 10
        ... chr1 50 55""", from_string=True)
        >>> fasta = pybedtools.example_filename('test.fa')
        >>> a = a.sequence(fi=fasta)
        >>> print open(a.seqfn).read()
        >chr1:1-10
        GATGAGTCT
        >chr1:50-55
        CCATC
        <BLANKLINE>

        '''
        kwargs['fi'] = fi

        if 'bed' not in kwargs:
            kwargs['bed'] = self.fn

        if 'fo' not in kwargs:
            tmp = self._tmp()
            kwargs['fo'] = tmp

        def check_sequence_stderr(x):
            if x.startswith('index file'):
                return True
            return False

        cmds, tmp, stdin = self.handle_kwargs(prog='fastaFromBed', **kwargs)
        _ = call_bedtools(cmds, tmp, stdin=stdin, \
                          check_stderr=check_sequence_stderr)
        self.seqfn = kwargs['fo']
        return self

    @_help('subtractBed')
    @_file_or_bedtool()
    @_returns_bedtool()
    @_log_to_history
    def subtract(self, b=None, **kwargs):
        """
        Subtracts from another BED file and returns a new BedTool object.

        Example usage:

            >>> a = pybedtools.example_bedtool('a.bed')
            >>> b = pybedtools.example_bedtool('b.bed')

            Do a "stranded" subtraction:

            >>> c = a.subtract(b, s=True)

            Require 50% of features in `a` to overlap:

            >>> c = a.subtract(b, f=0.5)
        """
        kwargs['b'] = b

        if 'a' not in kwargs:
            kwargs['a'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='subtractBed', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    @_help('slopBed')
    @_implicit('-i')
    @_returns_bedtool()
    @_log_to_history
    def slop(self, **kwargs):
        """
        Wraps slopBed, which adds bp to each feature.  Returns a new BedTool
        object.

        If *g* is a dictionary (for example, return values from
        pybedtools.chromsizes() ) it will be converted to a temp file for use
        with slopBed.  If it is a string, then it is assumed to be a filename.

        Alternatively, use *genome* to indicate a pybedtools-created genome.
        Example usage:

            >>> a = pybedtools.example_bedtool('a.bed')

            Increase the size of features by 100 bp in either direction.  Note
            that you need to specify either a dictionary of chromsizes or a
            filename containing chromsizes for the genome that your bed file
            corresponds to:

            >>> c = a.slop(g=pybedtools.chromsizes('hg19'), b=100)

            Grow features by 10 bp upstream and 500 bp downstream, using a
            genome file you already have constructed called 'hg19.genome'

            First, create the file:

            >>> fout = open('hg19.genome','w')
            >>> chromdict = pybedtools.get_chromsizes_from_ucsc('hg19')
            >>> for chrom, size in chromdict.items():
            ...     fout.write("%s\\t%s\\n" % (chrom, size[1]))
            >>> fout.close()

            Then use it:

            >>> c = a.slop(g='hg19.genome', l=10, r=500, s=True)

            Clean up afterwards:

            >>> os.unlink('hg19.genome')

        """
        if 'i' not in kwargs:
            kwargs['i'] = self.fn

        kwargs = self.check_genome(**kwargs)

        cmds, tmp, stdin = self.handle_kwargs(prog='slopBed', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    def check_genome(self, **kwargs):
        """
        Handles the different ways of specifying a genome in kwargs:

        g='genome.file' specifies a file directly
        genome='dm3' gets the file from genome registry
        self.chromsizes could be a dict.\
        """

        # If both g and genome are missing, assume self.chromsizes
        if ('g' not in kwargs) and ('genome' not in kwargs):
            if hasattr(self, 'chromsizes'):
                kwargs['g'] = self.chromsizes
            else:
                raise ValueError('No genome specified. Use the "g" or '
                                 '"genome" kwargs, or use the '
                                 '.set_chromsizes() method')

        # If both specified, rather than make an implicit decision, raise an
        # exception
        if 'g' in kwargs and 'genome' in kwargs:
            raise ValueError('Cannot specify both "g" and "genome"')

        # Something like genome='dm3' was specified
        if 'g' not in kwargs and 'genome' in kwargs:
            if isinstance(kwargs['genome'], dict):
                genome_dict = kwargs['genome']
            else:
                genome_dict = pybedtools.chromsizes(kwargs['genome'])
            genome_file = pybedtools.chromsizes_to_file(genome_dict)
            kwargs['g'] = genome_file
            del kwargs['genome']

        # By the time we get here, 'g' is specified.

        # If a dict was provided, convert to tempfile here
        if isinstance(kwargs['g'], dict):
            kwargs['g'] = pybedtools.chromsizes_to_file(kwargs['g'])

        if not os.path.exists(kwargs['g']):
            raise ValueError('Genome file "%s" does not exist')

        return kwargs

    @_help('mergeBed')
    @_implicit('-i')
    @_returns_bedtool()
    @_log_to_history
    def merge(self, **kwargs):
        """
        Merge overlapping features together. Returns a new BedTool object.

        Example usage:

            >>> a = pybedtools.example_bedtool('a.bed')

            Merge:

            >>> c = a.merge()

            Allow merging of features 500 bp apart:

            >>> c = a.merge(d=500)

            Report number of merged features:

            >>> c = a.merge(n=True)

            Report names of merged features:

            >>> c = a.merge(nms=True)

        """
        if 'i' not in kwargs:
            kwargs['i'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='mergeBed', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    @_help('closestBed')
    @_file_or_bedtool()
    @_implicit('-a')
    @_returns_bedtool()
    @_log_to_history
    def closest(self, b=None, **kwargs):
        """
        Return a new BedTool object containing closest features in *b*.  Note
        that the resulting file is no longer a valid BED format; use the
        special "_closest" methods to work with the resulting file.

        Example usage::

            a = BedTool('in.bed')

            # get the closest feature in 'other.bed' on the same strand
            b = a.closest('other.bed', s=True)

        """
        kwargs['b'] = b

        if ('abam' not in kwargs) and ('a' not in kwargs):
            kwargs['a'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='closestBed', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    @_help('windowBed')
    @_file_or_bedtool()
    @_implicit('-a')
    @_returns_bedtool()
    @_log_to_history
    def window(self, b=None, **kwargs):
        """
        Intersect with a window.

        Example usage::

            >>> a = pybedtools.example_bedtool('a.bed')
            >>> b = pybedtools.example_bedtool('b.bed')
            >>> print a.window(b, w=1000) #doctest: +NORMALIZE_WHITESPACE
            chr1	1	100	feature1	0	+	chr1	155	200	feature5	0	-
            chr1	1	100	feature1	0	+	chr1	800	901	feature6	0	+
            chr1	100	200	feature2	0	+	chr1	155	200	feature5	0	-
            chr1	100	200	feature2	0	+	chr1	800	901	feature6	0	+
            chr1	150	500	feature3	0	-	chr1	155	200	feature5	0	-
            chr1	150	500	feature3	0	-	chr1	800	901	feature6	0	+
            chr1	900	950	feature4	0	+	chr1	155	200	feature5	0	-
            chr1	900	950	feature4	0	+	chr1	800	901	feature6	0	+
            <BLANKLINE>
        """
        kwargs['b'] = b

        if ('abam' not in kwargs) and ('a' not in kwargs):
            kwargs['a'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='windowBed', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    @_help('shuffleBed')
    @_implicit('-i')
    @_log_to_history
    def shuffle(self, **kwargs):
        """
        Shuffle coordinates.

        Example usage:

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> seed = 1 # so this test always returns the same results
        >>> b = a.shuffle(genome='hg19', chrom=True, seed=seed)
        >>> print b #doctest: +NORMALIZE_WHITESPACE
        chr1	59535036	59535135	feature1	0	+
        chr1	99179023	99179123	feature2	0	+
        chr1	186189051	186189401	feature3	0	-
        chr1	219133189	219133239	feature4	0	+
        <BLANKLINE>

        """
        kwargs = self.check_genome(**kwargs)

        if 'i' not in kwargs:
            kwargs['i'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='shuffleBed', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    @_help('sortBed')
    @_implicit('-i')
    @_returns_bedtool()
    @_log_to_history
    def sort(self, **kwargs):
        """
        Note that chromosomes are sorted lexograpically, so chr12 will come
        before chr9.

        Example usage:

        >>> a = pybedtools.BedTool('''
        ... chr9 300 400
        ... chr1 100 200
        ... chr1 1 50
        ... chr12 1 100
        ... chr9 500 600
        ... ''', from_string=True)
        >>> print a.sort() #doctest: +NORMALIZE_WHITESPACE
        chr1	1	50
        chr1	100	200
        chr12	1	100
        chr9	300	400
        chr9	500	600
        <BLANKLINE>

        """
        if 'i' not in kwargs:
            kwargs['i'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='sortBed', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    @_help('annotateBed')
    @_implicit('-i')
    @_returns_bedtool()
    @_log_to_history
    def annotate(self, **kwargs):
        """
        Annotate this BedTool with a list of other files.
        Example usage:

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> b_fn = pybedtools.example_filename('b.bed')
        >>> print a.annotate(files=b_fn) #doctest: +NORMALIZE_WHITESPACE
        chr1	1	100	feature1	0	+	0.000000
        chr1	100	200	feature2	0	+	0.450000
        chr1	150	500	feature3	0	-	0.128571
        chr1	900	950	feature4	0	+	0.020000
        <BLANKLINE>
        """
        if 'i' not in kwargs:
            kwargs['i'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='annotateBed', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    @_help('flankBed')
    @_implicit('-i')
    @_returns_bedtool()
    @_log_to_history
    def flank(self, **kwargs):
        """
        Create flanking intervals on either side of input BED.

        Example usage:

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> print a.flank(genome='hg19', b=100) #doctest: +NORMALIZE_WHITESPACE
        chr1	0	1	feature1	0	+
        chr1	100	200	feature1	0	+
        chr1	0	100	feature2	0	+
        chr1	200	300	feature2	0	+
        chr1	50	150	feature3	0	-
        chr1	500	600	feature3	0	-
        chr1	800	900	feature4	0	+
        chr1	950	1050	feature4	0	+
        <BLANKLINE>

        """
        kwargs = self.check_genome(**kwargs)

        if 'i' not in kwargs:
            kwargs['i'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='flankBed', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    @_help('genomeCoverageBed')
    @_implicit('-i')
    @_returns_bedtool()
    @_log_to_history
    def genome_coverage(self, **kwargs):
        """
        Calculates coverage at each position in the genome.

        Use *bg=True* to have the resulting BedTool return valid BED-like
        features

        Example usage:

        >>> a = pybedtools.example_bedtool('x.bam')
        >>> b = a.genome_coverage(ibam=a.fn, genome='dm3', bg=True)
        >>> b.head(3) #doctest: +NORMALIZE_WHITESPACE
        chr2L	9329	9365	1
        chr2L	10212	10248	1
        chr2L	10255	10291	1

        """
        kwargs = self.check_genome(**kwargs)

        if 'i' not in kwargs:
            kwargs['i'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='genomeCoverageBed',
                                              **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    @_help('coverageBed')
    @_implicit('-a')
    @_returns_bedtool()
    @_log_to_history
    def coverage(self, b=None, **kwargs):
        """
        >>> a = pybedtools.example_bedtool('a.bed')
        >>> b = pybedtools.example_bedtool('b.bed')
        >>> c = a.coverage(b)
        >>> c.head(3) #doctest: +NORMALIZE_WHITESPACE
        chr1	155	200	feature5	0	-	2	45	45	1.0000000
        chr1	800	901	feature6	0	+	1	1	101	0.0099010
        """
        if ('abam' not in kwargs) and ('a' not in kwargs):
            kwargs['a'] = self.fn

        kwargs['b'] = b

        cmds, tmp, stdin = self.handle_kwargs(prog='coverageBed',
                                              **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

    @_help('maskFastaFromBed')
    @_implicit('-bed')
    @_log_to_history
    @_returns_bedtool()
    def mask_fasta(self, **kwargs):
        """
        Masks a fasta file at the positions in a BED file and saves result as
        *out*. This method returns None, and sets self.seqfn to *out*.

        >>> a = pybedtools.BedTool('chr1 100 110', from_string=True)
        >>> fasta_fn = pybedtools.example_filename('test.fa')
        >>> a = a.mask_fasta(fi=fasta_fn, fo='masked.fa.example')
        >>> b = a.slop(b=2, genome='hg19')
        >>> b = b.sequence(a.seqfn)
        >>> print b.print_sequence()
        >chr1:98-112
        TTNNNNNNNNNNAT
        <BLANKLINE>
        >>> os.unlink('masked.fa.example')
        >>> if os.path.exists('masked.fa.example.fai'):
        ...     os.unlink('masked.fa.example.fai')

        """
        if 'bed' not in kwargs:
            kwargs['bed'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='maskFastaFromBed',
                                              **kwargs)
        _ = call_bedtools(cmds, tmp, stdin=stdin)
        self.seqfn = kwargs['fo']
        return self

    def features(self):
        """
        Returns an iterator of :class:`feature` objects.
        """
        return iter(self)

    def count(self):
        """
        Number of features in BED file. Does the same thing as len(self), which
        actually just calls this method.

        Only counts the actual features.  Ignores any track lines, browser
        lines, lines starting with a "#", or blank lines.

        Example usage::

            a = BedTool('in.bed')
            a.count()
        """
        return sum(1 for _ in self)

    def print_sequence(self):
        """
        Print the sequence that was retrieved by the :meth:`BedTool.sequence`
        method.

        See usage example in :meth:`BedTool.sequence`.
        """
        if not hasattr(self, 'seqfn'):
            raise ValueError('Use .sequence(fasta) to get the sequence first')
        f = open(self.seqfn)
        s = f.read()
        f.close()
        return s

    def save_seqs(self, fn):
        """
        Save sequences of features in this BedTool object as a fasta file *fn*.

        In order to use this function, you need to have called
        the :meth:`BedTool.sequence()` method.

        A new BedTool object is returned which references the newly saved file.

        Example usage::

            a = BedTool('in.bed')

            # specify the filename of the genome in fasta format
            a.sequence('data/genomes/genome.fa')

            # use this method to save the seqs that correspond to the features
            # in "a"
            a.save_seqs('seqs.fa')
        """
        if not hasattr(self, 'seqfn'):
            raise ValueError('Use .sequence(fasta) to get the sequence first')
        fout = open(fn, 'w')
        fout.write(open(self.seqfn).read())
        fout.close()
        new_bedtool = BedTool(self.fn)
        new_bedtool.seqfn = fn
        return new_bedtool

    def randomstats(self, other, iterations, **kwargs):
        """
        Sends args and kwargs to :meth:`BedTool.randomintersection` and
        compiles results into a dictionary with useful stats.  Requires scipy
        and numpy.

        This is one possible way of assigning significance to overlaps between
        two files. See, for example:

            Negre N, Brown CD, Shah PK, Kheradpour P, Morrison CA, et al. 2010
            A Comprehensive Map of Insulator Elements for the Drosophila
            Genome. PLoS Genet 6(1): e1000814. doi:10.1371/journal.pgen.1000814

        Example usage:

        Make chromsizes a very small genome for this example:
        >>> chromsizes = {'chr1':(1,1000)}
        >>> a = pybedtools.example_bedtool('a.bed').set_chromsizes(chromsizes)
        >>> b = pybedtools.example_bedtool('b.bed')
        >>> results = a.randomstats(b, 100, debug=True)

        *results* is a dictionary that you can inspect.  The actual overlap:
        >>> print results['actual']
        3

        The median of all randomized overlaps:
        >>> print results['median randomized']
        2.0

        The percentile of the actual overlap in the distribution of randomized
        overlaps, which can be used to get an empirical p-value:
        >>> print results['percentile']
        90.0
        """
        if 'intersect_kwargs' not in kwargs:
            kwargs['intersect_kwargs'] = {'u': True}
        try:
            from scipy import stats
            import numpy as np
        except ImportError:
            raise ImportError("Need to install NumPy and SciPy for stats...")

        if isinstance(other, basestring):
            other = BedTool(other)
        else:
            assert isinstance(other, BedTool),\
                 'Either filename or another BedTool instance required'

        # Actual (unshuffled) counts.
        actual = len(self.intersect(other, **kwargs['intersect_kwargs']))

        # List of counts from randomly shuffled versions.
        # Length of counts == *iterations*.
        distribution = self.randomintersection(other, iterations=iterations,
                                               **kwargs)
        distribution = np.array(list(distribution))

        # Median of distribution
        med_count = np.median(distribution)

        n = float(len(distribution))

        frac_above = sum(distribution >= actual) / n
        frac_below = sum(distribution <= actual) / n

        normalized = actual / med_count

        lower_thresh = 2.5
        upper_thresh = 97.5
        lower = stats.scoreatpercentile(distribution, lower_thresh)
        upper = stats.scoreatpercentile(distribution, upper_thresh)

        actual_percentile = stats.percentileofscore(distribution, actual)
        d = {
        'iterations': iterations,
            'actual': actual,
            'file_a': self.fn,
            'file_b': other.fn,
             self.fn: len(self),
            other.fn: len(other),
              'self': len(self),
             'other': len(other),
        'frac randomized above actual': frac_above,
        'frac randomized below actual': frac_below,
        'median randomized': med_count,
        'normalized': normalized,
        'percentile': actual_percentile,
        'lower_%sth' % lower_thresh: lower,
        'upper_%sth' % upper_thresh: upper,
        }
        return d

    def randomintersection(self, other, iterations, intersect_kwargs=None,
                           shuffle_kwargs=None, debug=False):
        """
        Performs *iterations* shufflings of self, each time intersecting with
        *other*.

        Returns a generator of integers where each integer is the number of
        intersections of a shuffled file with *other*. This distribution can
        be used in downstream analysis for things like empirical p-values.

        *intersect_kwargs* and *shuffle_kwargs* are passed to self.intersect()
        and self.shuffle() respectively.  By default for intersect, u=True is
        specified -- but s=True might be a useful option for strand-specific
        work.

        Useful kwargs for *shuffle_kwargs* are chrom, excl, or incl.  If you
        use the "seed" kwarg, that seed will be used *each* time shuffleBed is
        called -- so all your randomization results will be identical for each
        iteration.  To get around this and to allow for tests, debug=True will
        set the seed to the iteration number.

        Example usage:

            >>> chromsizes = {'chr1':(0, 1000)}
            >>> a = pybedtools.example_bedtool('a.bed').set_chromsizes(chromsizes)
            >>> b = pybedtools.example_bedtool('b.bed')
            >>> results = a.randomintersection(b, 10, debug=True)
            >>> print list(results)
            [2, 2, 2, 0, 2, 3, 2, 1, 2, 3]

        """

        if shuffle_kwargs is None:
            shuffle_kwargs = {}
        if intersect_kwargs is None:
            intersect_kwargs = {'u': True}

        if not 'u' in intersect_kwargs:
            intersect_kwargs['u'] = True

        for i in range(iterations):
            if debug:
                shuffle_kwargs['seed'] = i
            tmp = self.shuffle(**shuffle_kwargs)
            tmp2 = tmp.intersect(other, **intersect_kwargs)
            yield len(tmp2)
            os.unlink(tmp.fn)
            os.unlink(tmp2.fn)
            del(tmp)
            del(tmp2)

    @_file_or_bedtool()
    @_returns_bedtool()
    def cat(self, other, postmerge=True, **kwargs):
        """
        Concatenates two BedTool objects (or an object and a file) and does an
        optional post-merge of the features.

        Use *postmerge=False* if you want to keep features separate.

        TODO:

            currently truncates at BED3 format!

        kwargs are sent to :meth:`BedTool.merge`.

        Example usage:

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> b = pybedtools.example_bedtool('b.bed')
        >>> print a.cat(b) #doctest: +NORMALIZE_WHITESPACE
        chr1	1	500
        chr1	800	950
        <BLANKLINE>

        """
        tmp = self._tmp()
        if isinstance(other, basestring):
            other = BedTool(other)
        else:
            assert isinstance(other, BedTool),\
                    'Either filename or another BedTool instance required'
        TMP = open(tmp, 'w')
        for f in self:
            TMP.write('%s\t%i\t%i\n' % (f.chrom, f.start, f.end))
        for f in other:
            TMP.write('%s\t%i\t%i\n' % (f.chrom, f.start, f.end))
        TMP.close()
        c = BedTool(tmp)
        if postmerge:
            d = c.merge(**kwargs)
            return d
        else:
            return c

    @_returns_bedtool()
    def saveas(self, fn, trackline=None):
        """
        Save BED file as a new file, adding the optional *trackline* to the
        beginning.

        Returns a new BedTool for the newly saved file.

        A newline is automatically added to the trackline if it does not
        already have one.

        Example usage:

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> b = a.saveas('other.bed')
        >>> b.fn
        'other.bed'
        >>> print b == a
        True

        >>> b = a.saveas('other.bed', trackline="name='test run' color=0,55,0")
        >>> open(b.fn).readline()
        "name='test run' color=0,55,0\\n"
        """
        fout = open(fn, 'w')
        if trackline is not None:
            fout.write(trackline.strip() + '\n')
        fout.write(str(self))
        fout.close()
        return BedTool(fn)

    @_returns_bedtool()
    def random_subset(self, n):
        '''
        Returns a new bedtools object containing a random subset of the
        features in this subset.

        Example usage:

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> b = a.random_subset(2)
        >>> len(b)
        2
        '''
        idxs = range(len(self))
        random.shuffle(idxs)
        idxs = idxs[:n]

        tmpfn = self._tmp()
        tmp = open(tmpfn, 'w')
        for i, f in enumerate(self):
            if i in idxs:
                tmp.write(str(f) + '\n')
        tmp.close()
        return BedTool(tmpfn)

    def total_coverage(self):
        """
        Returns the total number of bases covered by this BED file.  Does a
        self.merge() first to remove potentially multiple-counting bases.

        Example usage:

        >>> a = pybedtools.example_bedtool('a.bed')

        This does a merge() first, so this is what the total coverage is
        counting:

        >>> print a.merge() #doctest: +NORMALIZE_WHITESPACE
        chr1	1	500
        chr1	900	950
        <BLANKLINE>

        >>> print a.total_coverage()
        549
        """
        b = self.merge()
        total_bp = 0
        for feature in b.features():
            total_bp += len(feature)
        return total_bp

    @_returns_bedtool()
    def with_attrs(self, **kwargs):
        """
        Given arbitrary keyword arguments, turns the keys and values into
        attributes.  Useful for labeling BedTools at creation time.

        Example usage:

        >>> # add a "label" attribute to each BedTool
        >>> a = pybedtools.example_bedtool('a.bed')\
                                   .with_attrs(label='transcription factor 1')
        >>> b = pybedtools.example_bedtool('b.bed')\
                                   .with_attrs(label='transcription factor 2')
        >>> for i in [a, b]:
        ...     print i.count(), 'features for', i.label
        4 features for transcription factor 1
        2 features for transcription factor 2

        """
        for key, value in kwargs.items():
            setattr(self, key, value)
        return self
