import tempfile
import shutil
import subprocess
import operator
import os
import sys
import random
import string
from itertools import islice
import multiprocessing

from pybedtools.helpers import get_tempdir, _tags,\
    History, HistoryStep, call_bedtools, _flatten_list, \
    _check_sequence_stderr, isBAM, BEDToolsError, _call_randomintersect
from cbedtools import IntervalFile, IntervalIterator
import pybedtools

tempfile_prefix = 'pybedtools.'
tempfile_suffix = '.tmp'


_implicit_registry = {}
_other_registry = {}
_bam_registry = {}


def _wraps(prog=None, implicit=None, bam=None, other=None, uses_genome=False,
           make_tempfile_for=None, check_stderr=None, add_to_bedtool=None,
           nonbam=None, force_bam=False):
    """
    Do-it-all wrapper, to be used as a decorator.

    *prog* is the name of the BEDTools program that will be called.  The help
    for this program will also be added to the decorated method's docstring.

    *implicit* is the BEDTools program arg that should be filled in
    automatically.

    *bam* will disable the implicit substitution if *bam* is in the kwargs.
    This is typically 'abam' or 'ibam' if the program accepts BAM input.

    *other* is the BEDTools program arg that is passed in as the second input,
    if supported.  Within the semantics of BEDTools, the typical case will be
    that if implicit='a' then other='b'; if implicit='i' then other=None.

    *uses_genome*, if True, will check for 'g' and/or 'genome' args and
    retrieve the corresponding genome files as needed.

    *make_tempfile_for* is used for the sequence methods and indicates which
    kwarg should have a tempfile made for it if it's not provided ('fo' for the
    sequence methods)

    *check_stderr*, if not None, is a function that accepts a string (which
    will be anything written to stdout when calling the wrapped program).  This
    function should return True if the string is OK, and False if it should
    truly be considered an error.  This is needed for wrapping fastaFromBed,
    which will report to stderr that it's creating an index file.

    *add_to_bedtool* is used for sequence methods.  It is a dictionary mapping
    kwargs to attributes to be created in the resulting BedTool.  Typically it
    is {'fo':'seqfn'} which will add the resulting sequence name to the
    BedTool's .seqfn attribute. If *add_to_bedtool* is not None, then the
    returned BedTool will be *self* with the added attribute.

    *nonbam* is a kwarg that even if the input file was a BAM, the output will
    *not* be BAM format.  For example, the `-bed` arg for intersectBed will
    cause the output to be in BED format, not BAM.  If not None, this can be a
    string, a list of strings, or the special string "ALL", which means that
    the wrapped program will never return BAM output.

    *force_bam*, if True, will force the output to be BAM.  This is used for
    bedToBam.
    """
    not_implemented = False

    # Call the program with -h to get help, which prints to stderr.
    try:
        p = subprocess.Popen([prog, '-h'],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        help_str = p.communicate()[1]

        # underscores throw of ReStructuredText syntax of docstrings, so
        # replace 'em
        help_str = help_str.replace('_', '**')

        # indent
        help_str = help_str.split('\n')
        help_str = ['\t' + i for i in help_str]
        help_str = '\n'.join(help_str)

    # If the program can't be found, then we'll eventually replace the method
    # with a version that does nothing raise a NotImplementedError (plus a
    # helpful message).
    except OSError:
        help_str = '"%s" does not appear to be installed '\
                       'or on the path, so this method is '\
                       'disabled.  Please install a more recent '\
                       'version of BEDTools and re-import to '\
                       'use this method.' % prog
        not_implemented = True

    def decorator(func):
        """
        Accepts a function to be wrapped; discards the original and returns a
        new, rebuilt-from-scratch function based on the kwargs passed to
        _wraps().
        """
        # Register the implicit (as well as bam and other) args in the global
        # registry.  BedTool.handle_kwargs() will access these at runtime.  The
        # registry is keyed by program name (like intersectBed).
        _implicit_registry[prog] = implicit
        if other is not None:
            _other_registry[prog] = other
        if bam is not None:
            _bam_registry[prog] = bam

        # Here's where we replace an unable-to-be-found program's method with
        # one that only returns a NotImplementedError
        if not_implemented:
            def not_implemented_func(*args, **kwargs):
                raise NotImplementedError(help_str)
            return not_implemented_func

        def wrapped(self, *args, **kwargs):
            """
            A newly created function that will be returned by the _wraps()
            decorator
            """
            # Only one non-keyword argument is supported; this is then assumed
            # to be "other" (e.g., `-b` for intersectBed)
            if len(args) > 0:
                assert len(args) == 1
                kwargs[other] = args[0]

            # Should this function handle genome files?
            if uses_genome:
                kwargs = self.check_genome(**kwargs)

            # Add the implicit values to kwargs.  If the current BedTool is
            # BAM, it will automatically be passed to the appropriate
            # BAM-support arg (like `-abam`).  But this also allows the user to
            # explicitly specify the abam kwarg, which will override the
            # auto-substitution.
            # Note: here, `implicit` is something like "a"; `bam` is something
            # like "abam"
            if (implicit not in kwargs) and (bam not in kwargs):
                if not self._isbam:
                    kwargs[implicit] = self.fn
                else:
                    # It is a bam file.  If this program supports BAM as the
                    # first input, then we set it here
                    if bam is not None:
                        kwargs[bam] = self.fn

                    # Otherwise, BEDTools can't currently handle it, so raise
                    # an exception.
                    else:
                        raise BEDToolsError('"%s" currently can\'t handle BAM '
                                'input, please use bam_to_bed() first.' % prog)

            # For sequence methods, we may need to make a tempfile that will
            # hold the resulting sequence.  For example, fastaFromBed needs to
            # make a tempfile for 'fo' if no 'fo' was explicitly specified by
            # the user.
            if make_tempfile_for is not None:
                if make_tempfile_for not in kwargs:
                    kwargs[make_tempfile_for] = self._tmp()

            # At runtime, this will parse the kwargs, convert streams to
            # tempfiles if needed, and return all the goodies
            cmds, tmp, stdin = self.handle_kwargs(prog=prog, **kwargs)

            # Do the actual call
            stream = call_bedtools(cmds, tmp, stdin=stdin,
                                   check_stderr=check_stderr)
            result = BedTool(stream)

            # Post-hoc editing of the BedTool -- for example, this is used for
            # the sequence methods to add a `seqfn` attribute to the resulting
            # BedTool.
            if add_to_bedtool is not None:
                for kw, attr in add_to_bedtool.items():
                    value = kwargs[kw]
                    setattr(self, attr, value)
                    result = self

            # Decide whether the output is BAM format or not.
            result_is_bam = False

            # By default, if the current BedTool is BAM, then the result should
            # be, too.
            if self._isbam:
                result_is_bam = True

            # If nonbam is "ALL", then this method will never return BAM
            # output.
            if nonbam == 'ALL':
                result_is_bam = False

            # If any of the `nonbam` args are found in kwargs, then result is
            # not a BAM.  Side note: the _nonbam name mangling is necessary to
            # keep the nonbam arg passed into the original _wraps() decorator
            # in scope.
            if nonbam is not None and nonbam != 'ALL':
                if isinstance(nonbam, basestring):
                    _nonbam = [nonbam]
                else:
                    _nonbam = nonbam
                for i in _nonbam:
                    if i in kwargs:
                        result_is_bam = False
                        break

            if force_bam:
                result_is_bam = True

            result._isbam = result_is_bam
            return result

        # Now add the edited docstring (original Python doctring plus BEDTools
        # help) to the newly created method above
        if func.__doc__ is None:
            orig = ''
        else:
            orig = func.__doc__
        wrapped.__doc__ = orig + help_str

        # Add the original method's name to a new attribute so we can access it
        # when logging history
        wrapped._name = func.__name__

        return wrapped

    return decorator


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
             >>> assert 'a.bed' in example_files
             >>> a = pybedtools.example_bedtool('a.bed')

        """
        self._isbam = False
        self._bam_header = ""

        if not from_string:
            if isinstance(fn, BedTool):
                fn = fn.fn
            elif isinstance(fn, basestring):
                if not os.path.exists(fn):
                    raise ValueError('File "%s" does not exist' % fn)
                self._isbam = isBAM(fn)
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
        self._file_type = None
        self.history = History()

        if self._isbam and isinstance(self.fn, basestring):
            self._bam_header = ''.join(BAM(self.fn, header_only=True))
        else:
            self._bam_header = ""

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

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> a.field_count()
        6
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

    def introns(self, gene="gene", exon="exon"):
        """
        Given a BED12 or a GFF with exons, create a new `BedTool` with
        just introns.
        the output is a bed6 file with the score column (5) being one of
        'intron'/'utr5'/'utr3'
        """
        # iterate over all the features in the gene.
        s = self.sort()
        if self.file_type == "gff":
            exon_iter = BedTool((f for f in s if f[2] == gene)).saveas()
            gene_iter = BedTool((f for f in s if f[2] == exon)).saveas()

        elif self.file_type == "bed":
            if s.field_count() == 12:
                exon_iter = s.bed6().saveas()
                gene_iter = s.saveas()
            else:
                # TODO: bed6. groupby on name and find smallest start,
                # largest stop.
                exon_iter = s
                gene_iter = None
                raise NotImplementedError('.introns() only supported for bed12'
                                          'and GFF')

        else:
            raise NotImplementedError('.introns() only'
                            'supported for BED and GFF')

        fh = open(BedTool._tmp(), "w")

        # group on the name.
        exon_intervals = exon_iter.intervals
        for g in gene_iter:
            # search finds all, but we just want the ones that completely
            # overlap this gene.
            exons = [e for e in exon_intervals.search(g, same_strand=True)
                    if e.start >= g.start and e.end <= g.end]

            for i, exon in enumerate(exons):
                # 5' utr between gene start and first intron
                if i == 0 and exon.start > g.start:
                    utr = {"+": "utr5", "-": "utr3"}[g.strand]
                    print >>fh, "%s\t%i\t%i\t%s\t%s\t%s" % (g.chrom,
                                                   g.start,
                                                   exon.start,
                                                   g.name,
                                                   utr, g.strand)
                elif i == len(exons) - 1 and exon.end < g.end:
                    utr = {"+": "utr3", "-": "utr5"}[g.strand]
                    print >>fh, "%s\t%i\t%i\t%s\t%s\t%s" % (g.chrom,
                                                   exon.end,
                                                   g.end,
                                                   g.name,
                                                   utr, g.strand)
                elif i != len(exons) - 1:
                    istart = exon.end
                    iend = exons[i + 1].start
                    print >>fh, "%s\t%i\t%i\t%s\tintron\t%s" % \
                                              (g.chrom,
                                               istart, iend,
                                               g.name,
                                               g.strand)
        fh.close()
        return BedTool(fh.name)

    @property
    def file_type(self):
        """
        Return the type of the current file.  One of ('bed','vcf','gff', 'bam',
        'sam', 'empty').

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> a.file_type
        'bed'
        """
        if self._file_type is None:
            if not isinstance(self.fn, basestring):
                raise ValueError('Checking file_type not supported for '
                                 'non-file BedTools. Use .saveas() to '
                                 'save as a temp file first.')
            if self._isbam:
                self._file_type = 'bam'
            else:
                try:
                    self._file_type = IntervalFile(self.fn).file_type
                except StopIteration:
                    self._file_type = 'empty'
                except ValueError:
                    self._file_type = IntervalIterator(open(self.fn))\
                            .next().file_type
        return self._file_type

    def cut(self, indexes, stream=False):
        """
        Similar to unix `cut` except indexes are 0-based, must be a list
        and the columns are returned in the order requested.

        In addition, indexes can contain keys of the GFF/GTF attributes,
        in which case the values are returned. e.g. 'gene_name' will return the
        corresponding name from a GTF, or 'start' will return the start
        attribute of a BED Interval.

        See .with_column() if you need to do more complex operations.
        """
        if stream:
            return BedTool(([f[attr] for attr in indexes] for f in self))
        else:
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
        tmpfn = tempfile.NamedTemporaryFile(prefix=tempfile_prefix,
                                            suffix=tempfile_suffix,
                                            delete=False)
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
            if self._isbam:
                # Note: BAM class takes filename or stream, so self.fn is OK
                # here
                return IntervalIterator(BAM(self.fn))

            # TODO: Sort of a hack, cause we can't use IntervalFile as a SAM
            # iterator [yet]
            elif self.file_type == 'sam':
                return IntervalIterator(open(self.fn))

            # Easy case: BED/GFF/VCF, as a file
            else:
                return IntervalFile(self.fn)

        # Open file, like subprocess.PIPE.
        if isinstance(self.fn, file):
            if self._isbam:
                return IntervalIterator(BAM(self.fn))
            else:
                # Note: even if this is a SAM, the filetype handling eventually
                # gets passed to create_interval_from_fields.
                return IntervalIterator(self.fn)

        # Otherwise assume fn is already an iterable
        else:
            return self.fn

    @property
    def intervals(self):
        return IntervalFile(self.fn)

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
        if isinstance(self.fn, basestring) and not self._isbam:
            return open(self.fn).read()

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
        if str(self) == str(other):
            return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __getitem__(self, key):
        if isinstance(key, slice):
            return islice(self, key.start, key.stop, key.step)
        elif isinstance(key, int):
            return list(islice(self, key, key + 1))[0]
        else:
            raise ValueError('Only slices or integers allowed for indexing '
                             'into a BedTool')

    def __add__(self, other):
        return self.intersect(other, u=True)

    def __sub__(self, other):
        return self.intersect(other, v=True)

    def head(self, n=10):
        """
        Prints the first *n* lines

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> a.head(2) #doctest: +NORMALIZE_WHITESPACE
        chr1	1	100	feature1	0	+
        chr1	100	200	feature2	0	+
        <BLANKLINE>

        """
        if not isinstance(self.fn, basestring):
            raise NotImplementedError('head() not supported for non file-based'
                    'BedTools')
        for i, line in enumerate(iter(self)):
            if i == (n):
                break
            print line

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

    def _collapse(self, iterable, fn=None, trackline=None):
        """
        Collapses an iterable into file *fn* (or a new tempfile if *fn* is
        None).

        Returns the newly created filename.
        """
        if fn is None:
            fn = self._tmp()

        fout = open(fn, 'w')

        if trackline:
            fout.write(trackline.strip() + '\n')

        for i in iterable:
            fout.write(str(i) + '\n')
        fout.close()
        return fn

    def handle_kwargs(self, prog, debug=False, **kwargs):
        """
        Handle most cases of BEDTool program calls, but leave the specifics
        up to individual methods.

        *prog* is a BEDTools program name, e.g., 'intersectBed'.

        If *debug*, prints kwargs recieved and cmds returned to stderr.

        *kwargs* are passed directly from the calling method (like
        self.intersect).

        This method figures out, given how this BedTool was constructed, what
        to send to BEDTools programs -- for example, an open file to stdin with
        the `-` argument, or a filename with the `-a` argument.
        """
        if debug:
            sys.stderr.write('kwargs: ' + str(kwargs) + '\n')

        # If you pass in a list, how should it be converted to a BedTools arg?
        default_list_delimiter = ' '
        list_delimiters = {'annotateBed': ' ',
                               'overlap': ',',
                               'groupBy': ','}

        stdin = None

        # -----------------------------------------------------------------
        # Decide how to send instream1 to BEDTools.  If there's no implicit
        # instream1 arg, then do nothing.
        #
        try:
            # e.g., 'a' for intersectBed
            if self._isbam:
                inarg1 = _bam_registry[prog]
            else:
                inarg1 = _implicit_registry[prog]

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
            inarg2 = _other_registry[prog]

            # e.g., another BedTool
            instream2 = kwargs[inarg2]

            # Get stream if BedTool
            if isinstance(instream2, BedTool):
                instream2 = instream2.fn

            # Filename
            if isinstance(instream2, basestring):
                kwargs[inarg2] = instream2

            # Otherwise we need to collapse it in order to send to BEDTools
            # programs
            else:
                kwargs[inarg2] = self._collapse(instream2)

        except KeyError:
            pass

        # If stream not specified, then a tempfile will be created
        if kwargs.pop('stream', None):
                tmp = None
        else:
            output = kwargs.pop('output', None)
            if output:
                tmp = output
            else:
                tmp = self._tmp()

        additional_args = kwargs.pop('additional_args', None)

        # Parse the kwargs into BEDTools-ready args
        cmds = [prog]
        for key, value in kwargs.items():
            if isinstance(value, bool):
                if value:
                    cmds.append('-' + key)
                else:
                    continue
            elif isinstance(value, list) or isinstance(value, tuple):
                value = map(str, value)
                try:
                    delim = list_delimiters[prog]
                except KeyError:
                    delim = default_list_delimiter

                if delim == ' ':
                    cmds.append('-' + key)
                    cmds.extend(value)

                # make comma-separated list if that's what's needed
                else:
                    cmds.append('-' + key)
                    cmds.append(delim.join(value))

            else:
                cmds.append('-' + key)
                cmds.append(str(value))

        if additional_args:
            cmds.append(additional_args)

        if debug:
            sys.stderr.write('cmds: ' + str(cmds) + '\n')
            sys.stderr.write('tmp: ' + repr(tmp) + '\n')
            sys.stderr.write('stdin: ' + repr(stdin) + '\n')
        return cmds, tmp, stdin

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

        # If it's a file-based BedTool -- which is likely, if we're trying to
        # remove invalid features -- then we need to parse it line by line.
        if isinstance(self.fn, basestring):
            i = IntervalIterator(open(self.fn))
        else:
            i = IntervalIterator(self.fn)

        def _generator():
            while True:
                try:
                    feature = i.next()
                    if feature.start < feature.stop:
                        yield feature
                except pybedtools.MalformedBedLineError:
                    continue
                except OverflowError:
                    # This can happen if coords are negative
                    continue
                except StopIteration:
                    break
        return BedTool(_generator())

    def all_hits(self, interval, same_strand=False, overlap=0.0):
        """
        Calls the `all_hits` method of an IntervalFile to return all intervals
        in this current BedTool that overlap `interval`.

        Notes:
                If this current BedTool is generator-based, it will be
                converted into a file first.

                If this current BedTool refers to a BAM file, it will be
                converted to a BED file first using default arguments.  If you
                don't want this to happen, please convert to BED first before
                using this method.
        """
        fn = self.fn
        if not isinstance(fn, basestring):
            fn = self.saveas().fn
        if self._isbam:
            fn = self.bam_to_bed().fn
        interval_file = pybedtools.IntervalFile(fn)
        return interval_file.all_hits(interval, same_strand, overlap)

    @_log_to_history
    @_wraps(prog='bed12ToBed6', implicit='i', bam=None, other=None)
    def bed6(self, **kwargs):
        """
        convert a BED12 to a BED6 file
        """
        pass

    @_log_to_history
    @_wraps(prog='bamToBed', implicit='i', other=None, nonbam='ALL', bam='i')
    def bam_to_bed(self, **kwargs):
        """
        Convert BAM to BED.
        """

    @_wraps(prog='bedToBam', implicit='i', uses_genome=True, force_bam=True)
    def _bed_to_bam(self):
        """
        Wraps bedToBam and is called internally for BED/GFF/VCF files by
        self.to_bam (which needs to do something different for SAM files...)
        """

    @_log_to_history
    def to_bam(self, **kwargs):
        """
        If self.fn is in BED/VCF/GFF format, call BEDTools' bedToBam.  If
        self.fn is in SAM format, then create a header out of the genome file
        and then convert using `samtools`.
        """
        if self.file_type in ('bed', 'gff', 'vcf'):
            return self._bed_to_bam(**kwargs)
        if self.file_type == 'sam':

            # construct a genome out of whatever kwargs were passed in
            kwargs = self.check_genome(**kwargs)

            cmds = [os.path.join(pybedtools._samtools_path, 'samtools'),
                    'view',
                    '-S',
                    '-b',
                    '-t', kwargs['g'],
                    '-']
            tmp = self._tmp()
            try:
                p = subprocess.Popen(cmds,
                                     stdout=open(tmp, 'w'),
                                     stderr=subprocess.PIPE,
                                     stdin=subprocess.PIPE,
                                     bufsize=1)
            except OSError:
                raise OSError('SAMtools (http://samtools.sourceforge.net/) '
                              'needs to be installed for BAM support')
            for line in self:
                p.stdin.write(str(line) + '\n')
            stdout, stderr = p.communicate()
            new_bedtool = BedTool(tmp)
            new_bedtool._isbam = True
            return new_bedtool

    @_log_to_history
    @_wraps(prog='intersectBed', implicit='a', other='b', bam='abam',
            nonbam='bed')
    def intersect(self):
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

    @_log_to_history
    @_wraps(prog='fastaFromBed', implicit='bed', bam=None, other='fi',
            make_tempfile_for='fo', check_stderr=_check_sequence_stderr,
            add_to_bedtool={'fo': 'seqfn'})
    def sequence(self):
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

    @_log_to_history
    @_wraps(prog='nucBed', implicit='bed', other='fi')
    def nucleotide_content(self):
        """
        Profiles nucleotide content.  The returned BED file contains extra
        information about the nucleotide content
        """

    @_log_to_history
    @_wraps(prog='multiBamCov', implicit='bed')
    def multi_bam_coverage(self):
        """
        Pass a list of sorted and indexed BAM files as `bams`
        """

    @_log_to_history
    @_wraps(prog='subtractBed', implicit='a', other='b', bam=None)
    def subtract(self):
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

    @_log_to_history
    @_wraps(prog='slopBed', implicit='i', other=None, bam=None,
            uses_genome=True)
    def slop(self):
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

    @_log_to_history
    @_wraps(prog='mergeBed', implicit='i', other=None, bam=None)
    def merge(self):
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

    @_log_to_history
    @_wraps(prog='closestBed', implicit='a', other='b', bam=None)
    def closest(self):
        """
        Return a new BedTool object containing closest features in *b*.  Note
        that the resulting file is no longer a valid BED format; use the
        special "_closest" methods to work with the resulting file.

        Example usage::

            a = BedTool('in.bed')

            # get the closest feature in 'other.bed' on the same strand
            b = a.closest('other.bed', s=True)

        """

    @_log_to_history
    @_wraps(prog='windowBed', implicit='a', other='b', bam='abam',
            nonbam='bed')
    def window(self):
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

    @_log_to_history
    @_wraps(prog='shuffleBed', implicit='i', other=None, bam=None,
            uses_genome=True)
    def shuffle(self):
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

    @_log_to_history
    @_wraps(prog='sortBed', implicit='i')
    def sort(self):
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

    @_log_to_history
    @_wraps(prog='annotateBed', implicit='i')
    def annotate(self):
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

    @_log_to_history
    @_wraps(prog='flankBed', implicit='i', uses_genome=True)
    def flank(self):
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

    @_log_to_history
    @_wraps(prog='genomeCoverageBed', implicit='i', bam='ibam',
            uses_genome=True, nonbam='ALL')
    def genome_coverage(self):
        """
        Calculates coverage at each position in the genome.

        Use *bg=True* to have the resulting BedTool return valid BED-like
        features

        Example usage:

        >>> a = pybedtools.example_bedtool('x.bam')
        >>> b = a.genome_coverage(genome='dm3', bg=True)
        >>> b.head(3) #doctest: +NORMALIZE_WHITESPACE
        chr2L	9329	9365	1
        chr2L	10212	10248	1
        chr2L	10255	10291	1
        """

    @_log_to_history
    @_wraps(prog='coverageBed', implicit='a', other='b', bam='abam',
            nonbam='ALL')
    def coverage(self):
        """
        >>> a = pybedtools.example_bedtool('a.bed')
        >>> b = pybedtools.example_bedtool('b.bed')
        >>> c = a.coverage(b)
        >>> c.head(3) #doctest: +NORMALIZE_WHITESPACE
        chr1	155	200	feature5	0	-	2	45	45	1.0000000
        chr1	800	901	feature6	0	+	1	1	101	0.0099010
        """

    @_log_to_history
    @_wraps(prog='maskFastaFromBed', implicit='bed', other='fi',
            make_tempfile_for='fo', add_to_bedtool={'fo': 'seqfn'},
            check_stderr=_check_sequence_stderr)
    def mask_fasta(self):
        """
        Masks a fasta file at the positions in a BED file and saves result as
        'out' and stores the filename in seqfn.

        >>> a = pybedtools.BedTool('chr1 100 110', from_string=True)
        >>> fasta_fn = pybedtools.example_filename('test.fa')
        >>> a = a.mask_fasta(fi=fasta_fn, fo='masked.fa.example')
        >>> b = a.slop(b=2, genome='hg19')
        >>> b = b.sequence(fi=a.seqfn)
        >>> print open(b.seqfn).read()
        >chr1:98-112
        TTNNNNNNNNNNAT
        <BLANKLINE>
        >>> os.unlink('masked.fa.example')
        >>> if os.path.exists('masked.fa.example.fai'):
        ...     os.unlink('masked.fa.example.fai')
        """

    @_log_to_history
    @_wraps(prog='complementBed', implicit='i', uses_genome=True)
    def complement(self):
        """
        >>> a = pybedtools.example_bedtool('a.bed')
        >>> a.complement(genome='hg19').head(5) #doctest: +NORMALIZE_WHITESPACE
        chr1	0	1
        chr1	500	900
        chr1	950	249250621
        chr10	0	135534747
        chr11	0	135006516
        """

    @_log_to_history
    @_wraps(prog='overlap', implicit='i')
    def overlap(self):
        """
        >>> a = pybedtools.example_bedtool('a.bed')
        >>> b = pybedtools.example_bedtool('b.bed')
        >>> c = a.window(b, w=10).overlap(cols=[2,3,8,9])
        >>> print c #doctest: +NORMALIZE_WHITESPACE
        chr1	100	200	feature2	0	+	chr1	155	200	feature5	0	-	45
        chr1	150	500	feature3	0	-	chr1	155	200	feature5	0	-	45
        chr1	900	950	feature4	0	+	chr1	800	901	feature6	0	+	1
        <BLANKLINE>
        """

    # TODO: needs test files and doctests written
    @_log_to_history
    @_wraps(prog='pairToBed', implicit='a', other='b', bam='abam',
            nonbam='bedpe')
    def pair_to_bed(self):
        """
        """

    # TODO: needs test files and doctests written
    @_log_to_history
    @_wraps(prog='pairToPair', implicit='a', other='b')
    def pair_to_pair(self):
        """
        """

    @_log_to_history
    @_wraps(prog='groupBy', implicit='i')
    def groupby(self):
        """
        >>> a = pybedtools.example_bedtool('gdc.gff')
        >>> b = pybedtools.example_bedtool('gdc.bed')
        >>> c = a.intersect(b, c=True)
        >>> d = c.groupby(g=[1, 4, 5], c=10, ops=['sum'])
        >>> print d #doctest: +NORMALIZE_WHITESPACE
        chr2L	41	70	0
        chr2L	71	130	2
        chr2L	131	170	4
        chr2L	171	200	0
        chr2L	201	220	1
        chr2L	41	130	2
        chr2L	171	220	1
        chr2L	41	220	7
        chr2L	161	230	6
        chr2L	41	220	7
        <BLANKLINE>

        """

    @_log_to_history
    @_wraps(prog='tagBam', implicit='i', bam='i')
    def tag_bam(self):
        """
        `files` and `labels` should lists of equal length.

        """

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

        Example usage:

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> a.count()
        4
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

        Example usage:

        >>> a = pybedtools.BedTool('''
        ... chr1 1 10
        ... chr1 50 55''', from_string=True)
        >>> fasta = pybedtools.example_filename('test.fa')
        >>> a = a.sequence(fi=fasta)
        >>> print open(a.seqfn).read()
        >chr1:1-10
        GATGAGTCT
        >chr1:50-55
        CCATC
        <BLANKLINE>
        >>> b = a.save_seqs('example.fa')
        >>> assert open(b.fn).read() == open(a.fn).read()
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
        >>> try:
        ...     results = a.randomstats(b, 100, debug=True)
        ... except ImportError:
        ...     # allow doctests to pass if SciPy not installed
        ...     pass

        *results* is a dictionary that you can inspect.

        (Note that the following examples are not run as part of the doctests
        to avoid forcing users to install SciPy just to pass tests)

        The actual overlap::

            print results['actual']
            3

        The median of all randomized overlaps::

            print results['median randomized']
            2.0

        The percentile of the actual overlap in the distribution of randomized
        overlaps, which can be used to get an empirical p-value::

            print results['percentile']
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
                           shuffle_kwargs=None, debug=False,
                           report_iterations=False, processes=None,
                           _orig_processes=None):
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
        set the seed to the iteration number.  You may also break up the
        intersections across multiple processes with *processes* > 1.

        Example usage:

            >>> chromsizes = {'chr1':(0, 1000)}
            >>> a = pybedtools.example_bedtool('a.bed')
            >>> a = a.set_chromsizes(chromsizes)
            >>> b = pybedtools.example_bedtool('b.bed')
            >>> results = a.randomintersection(b, 10, debug=True)
            >>> print list(results)
            [2, 2, 2, 0, 2, 3, 2, 1, 2, 3]

        """
        if processes is not None:
            p = multiprocessing.Pool(processes)
            iterations_each = [iterations / processes] * processes
            iterations_each[-1] += iterations % processes
            results = [p.apply_async(_call_randomintersect, (self, other, it),
                              dict(intersect_kwargs=intersect_kwargs,
                                   shuffle_kwargs=shuffle_kwargs,
                                   debug=debug,
                                   report_iterations=report_iterations,
                                   _orig_processes=processes))
                 for it in iterations_each]
            for r in results:
                for value in r.get():
                    yield value
            raise StopIteration

        if shuffle_kwargs is None:
            shuffle_kwargs = {}
        if intersect_kwargs is None:
            intersect_kwargs = {'u': True}

        if not 'u' in intersect_kwargs:
            intersect_kwargs['u'] = True

        resort = intersect_kwargs.get('sorted', False)

        for i in range(iterations):
            if debug:
                shuffle_kwargs['seed'] = i
            if report_iterations:
                if _orig_processes > 1:
                    msg = '\rapprox (total across %s processes): %s' \
                            % (_orig_processes, i * _orig_processes)
                else:
                    msg = '\r%s' % i
                sys.stderr.write(msg)
                sys.stderr.flush()

            # Re-sort if sorted=True in kwargs
            if resort:
                tmp0 = self.shuffle(stream=True, **shuffle_kwargs)
                tmp = tmp0.sort(stream=True)
            else:
                tmp = self.shuffle(stream=True, **shuffle_kwargs)

            tmp2 = tmp.intersect(other, stream=True, **intersect_kwargs)

            yield len(tmp2)

            # Close the open stdouts from subprocess.Popen calls.  Note: doing
            # this in self.__del__ doesn't fix the open file limit bug; it
            # needs to be done here.
            if resort:
                tmp0.fn.close()
            tmp.fn.close()
            tmp2.fn.close()
            del(tmp)
            del(tmp2)

    @_log_to_history
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

    @_log_to_history
    def saveas(self, fn=None, trackline=None):
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
        if fn is None:
            fn = self._tmp()

        fn = self._collapse(self, fn=fn, trackline=trackline)
        return BedTool(fn)

    @_log_to_history
    def moveto(self, fn=None):
        """
        Move BED file to new filename, `fn`.

        Returns a new BedTool for the new file.

        Example usage:

        >>> # make a copy so we don't mess up the example file
        >>> a = pybedtools.example_bedtool('a.bed').saveas()
        >>> a_contents = str(a)
        >>> b = a.moveto('other.bed')
        >>> b.fn
        'other.bed'
        >>> b == a_contents
        True
        """
        if not isinstance(self.fn, basestring):
            fn = self._collapse(self, fn=fn)
        else:
            shutil.move(self.fn, fn)
        return BedTool(fn)

    @_log_to_history
    def random_subset(self, n, seed=None):
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
        if seed is not None:
            random.seed(seed)
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

    @_log_to_history
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

    def as_intervalfile(self):
        """
        Returns an IntervalFile of this BedTool, which provides a low-level
        interface
        """
        if not isinstance(self.fn, basestring):
            fn = self._collapse(self.fn)
        else:
            fn = self.fn
        return IntervalFile(fn)


class BAM(object):
    def __init__(self, stream, header_only=False):
        """
        Wraps samtools to iterate over a BAM, yielding lines.
        """
        self.stream = stream
        self.header_only = header_only
        test_cmds = ['samtools', 'view']
        try:
            p = subprocess.Popen(test_cmds,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()
        except OSError:
            raise OSError('SAMtools (http://samtools.sourceforge.net/) '
                          'needs to be installed for BAM support')

        if isinstance(self.stream, basestring):
            self.cmds = [os.path.join(pybedtools._samtools_path, 'samtools'),
                         'view', stream]
            if header_only:
                self.cmds.append('-h')
            self.p = subprocess.Popen(self.cmds,
                                      stdout=subprocess.PIPE,
                                      bufsize=1)
        else:
            # Streaming . . .
            self.cmds = [os.path.join(pybedtools._samtools_path, 'samtools'),
                         'view', '-']
            if header_only:
                self.cmds.append('-h')
            self.p = subprocess.Popen(self.cmds,
                                      stdout=subprocess.PIPE,
                                      stdin=subprocess.PIPE,
                                      bufsize=0)
            # Can't iterate (for i in stream) cause we're dealing with a binary
            # BAM file here.  So read the whole thing in at once.
            self.p.stdin.write(stream.read())

    def __iter__(self):
        return self

    def next(self):
        line = self.p.stdout.next()

        # If we only want the header, then short-circuit once we're out of
        # header lines
        if self.header_only:
            if line[0] != '@':
                raise StopIteration

        return line
