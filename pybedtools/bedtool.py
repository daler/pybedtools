import sys
import os
import tempfile
import subprocess
import random
import itertools
import glob
from math import floor, ceil

from features import bedfeature
import genome_registry

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

def cleanup(verbose=False,remove_all=False):
    """
    Deletes all temporary files in the *bedtool.TEMPFILES* class
    variable.
    
    If *verbose*, reports what it's doing

    If *remove_all*, then ALL files matching "pybedtools.*.tmp" in the temp dir
    will be deleted.
    """
    for fn in bedtool.TEMPFILES:
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

class bedtool(object):
    """
    Wrapper around Aaron Quinlans ``BEDtools`` suite of programs
    (https://github.com/arq5x/bedtools); also contains many useful
    methods for more detailed work with BED files.

    Typical usage is to point to an existing file::

        a = bedtool('a.bed')

    But you can also create one from scratch from a string::

        s = '''
            chrX  1  100
            chrX 25  800
            '''
        a = bedtool(s,from_string=True).saveas('a.bed')

    Or use examples that come with pybedtools::
        
        example_beds = pybedtools.list_example_beds()
        a = pybedtools.example_bedtool('a.bed')


    """
    TEMPFILES = []
    def __init__(self,fn,genome=None, from_string=False):
        """
        *fn* is a BED format file, or alternatively another bedtool instance.

        *genome* is an optional genome assembly ('dm3', 'hg18', etc) or a
        dictionary of chrom:(start,stop) integers to consider as the genome
        space.  This is used for randomizations and coverage.

        If *from_string* is True, then treat all spaces as TABs and write
        to tempfile, treating *fn* as the contents of the bed file.  
        This strips empty lines.
        """
        if not from_string:
            if isinstance(fn, bedtool):
                fn = fn.fn
            if not os.path.exists(fn):
                raise ValueError, 'File "%s" does not exist' % fn
        else:
            bed_contents = fn
            fn = self._tmp()
            fout = open(fn,'w')
            for line in bed_contents.splitlines():
                if len(line) == 0:
                    continue
                line = '\t'.join(line.split())+'\n'
                fout.write(line)
            fout.close()

        self.fn = fn
        self._hascounts = False
        if genome is None:
            self.genome = None

        elif isinstance(genome,basestring):
            try:
                self.genome = getattr(genome_registry,genome)
            except KeyError:
                raise ValueError, 'Genome %s not registered' % genome
        else:
            self.genome = genome

    def _tmp(self):
        '''
        Makes a tempfile and registers it the the bedtool.TEMPFILES class variable.
        Adds a "pybedtools." prefix and ".tmp" extension for easy deletion if
        you forget to call pybedtools.cleanup().
        '''

        tmpfn = tempfile.NamedTemporaryFile(prefix='pybedtools.',suffix='.tmp',delete=False)
        tmpfn = tmpfn.name
        bedtool.TEMPFILES.append(tmpfn)
        return tmpfn

    def __iterator(self):
        '''Iterator that returns lines from BED file'''
        f = open(self.fn)
        for line in f:
            if line.startswith('browser'):
                continue
            if line.startswith('track'):
                continue
            if line.startswith('#'):
                continue
            if len(line.strip()) == 0:
                continue
            yield line
        f.close()
    
    def __iter__(self):
        return self.__iterator()

    def __repr__(self):
        return '<bedtool (%s)>'%self.fn

    def __str__(self):
        f = open(self.fn)
        s = f.read()
        f.close()
        return s

    def __len__(self):
        return self.count()

    @_file_or_bedtool()
    def __add__(self,other):
        return self.intersect(other,u=True)

    @_file_or_bedtool()
    def __sub__(self,other):
        return self.intersect(other, v=True)
        
    @property
    def lines(self):
        return open(self.fn)
   
    def head(self,n=10):
        """
        Prints the first *n* lines
        """
        for i,line in enumerate(self.lines):
            if i == (n):
                break
            print line,

    @_help('intersectBed')
    @_file_or_bedtool()
    @_implicit('-a')
    @_returns_bedtool()
    def intersect(self, b=None, **kwargs):
        """
        Intersect with another BED file. If you want to use BAM as input, you
        need to specify *abam='filename.bam'*.  Returns a new bedtool object.
       
        Example usage::

            # create new bedtool object
            a = bedtool('in.bed')

            # get overlaps with "other.bed"
            overlaps = a.intersect('other.bed')

            # use v=True to get the inverse, or those unique to in.bed
            unique_to_a = a.intersect('other.bed', v=True)
            
            # features unique to "other.bed"
            unique_to_other = bedtool('other.bed').intersect(a, v=True)

        """
        other = b
        if (type(other) is str) or (type(other) is unicode):
            kwargs['b'] = other
        else: 
            assert isinstance(other,bedtool), 'Either filename or another bedtool instance required'
            kwargs['b'] = other.fn
            
        tmp = self._tmp()
        cmds = ['intersectBed',]
        if 'abam' not in kwargs:
            kwargs['a'] = self.fn
        cmds.extend(self.parse_kwargs(**kwargs))
        cmds.extend(['>',tmp])
        os.system(' '.join(cmds))
        other = bedtool(tmp)
        if 'c' in kwargs:
            other._hascounts = True
        return other

    @_help('fastaFromBed')
    @_returns_bedtool()
    def sequence(self, **kwargs):
        '''
        Wraps ``fastaFromBed``.  *fi* is passed in by the user; *bed* is
        automatically passed in as the bedfile of this object; *fo* by default
        is a temp file.  Use save_seqs() to save as a file.
        
        Example usage::

            a = bedtool('in.bed')
            a.sequence(fi='genome.fa')
            a.print_sequence()
        '''
        tmp = self._tmp()
        kwargs['bed'] = self.fn
        kwargs['fo'] = tmp
        cmds = ['fastaFromBed']
        cmds.extend(self.parse_kwargs(**kwargs))
        os.system(' '.join(cmds))
        self.seqfn = tmp
        return self

    @_help('subtractBed')
    @_file_or_bedtool()
    @_returns_bedtool()
    def subtract(self, other, **kwargs):
        """
        Subtracts from another BED file and returns a new bedtool object.

        Example usage::

            a = bedtool('in.bed')

            # do a "stranded" subtraction
            b = a.subtract('other.bed',s=True)

            # Require 50% of features in a to overlap
            c = a.subtract('other.bed', s=0.5)

        """
        kwargs['a'] = self.fn
        if (type(other) is str) or (type(other) is unicode):
            kwargs['b'] = other
        else:
            assert isinstance(other,bedtool), 'Either filename or another bedtool instance required'
            kwargs['b'] = other.fn
        tmp = self._tmp()
        cmds = ['subtractBed',]
        cmds.extend(self.parse_kwargs(**kwargs))
        cmds.extend(['>',tmp])
        os.system(' '.join(cmds))
        return bedtool(tmp)

    @_help('slopBed')
    @_implicit('-i')
    @_returns_bedtool()
    def slop(self, genome=None, **kwargs):
        """
        Wraps slopBed, which adds bp to each feature.  Returns a new bedtool
        object.
        
        If *genome* is specified with a genome name, the genome file will be
        automatically retrieved from UCSC Genome Browser.

        Example usage::

            a = bedtool('in.bed')

            # increase the size of features by 100 bp in either direction
            b = a.slop(genome='dm3', b=100)

            # grow features by 10 bp upstream and 500 bp downstream,
            # using a genome file you already have constructed called
            # dm3.genome.
            c = a.slop(g='dm3.genome', l=10, r=500, s=True) 
        """
        if genome is not None:
            genome_fn = self.get_genome(genome)
            kwargs['g'] = genome_fn
        kwargs['i'] = self.fn
        tmp = self._tmp()
        cmds = ['slopBed',]
        cmds.extend(self.parse_kwargs(**kwargs))
        cmds.extend(['>',tmp])
        os.system(' '.join(cmds))
        return bedtool(tmp)

    @_help('mergeBed')
    @_implicit('-i')
    @_returns_bedtool()
    def merge(self, **kwargs):
        """
        Merge overlapping features together. Returns a new bedtool object.

        Example usage::

            a = bedtool('in.bed')

            # allow merging of features 100 bp apart
            b = a.merge(d=100)

        """
        tmp = self._tmp()
        cmds = ['mergeBed',]
        kwargs['i'] = self.fn
        cmds.extend(self.parse_kwargs(**kwargs))
        cmds.extend(['>',tmp])
        os.system(' '.join(cmds))
        return bedtool(tmp)

    @_help('closestBed')
    @_file_or_bedtool()
    @_implicit('-a')
    @_returns_bedtool()
    def closest(self, other, **kwargs):
        """
        Return a new bedtool object containing closest features in *other*.  Note
        that the resulting file is no longer a valid BED format; use the
        special "_closest" methods to work with the resulting file.

        Example usage::

            a = bedtool('in.bed')

            # get the closest feature in 'other.bed' on the same strand
            b = a.closest('other.bed', s=True)

        """
        tmp = self._tmp()
        cmds = ['closestBed',]
        kwargs['a'] = self.fn
        if (type(other) is str) or (type(other) is unicode):
            kwargs['b'] = other
        else:
            assert isinstance(other,bedtool), 'Either filename or another bedtool instance required'
            kwargs['b'] = other.fn
        cmds.extend(self.parse_kwargs(**kwargs))
        cmds.extend(['>',tmp])
        os.system(' '.join(cmds))
        newbedtool = bedtool(tmp)
        newbedtool.closest_output = True
        return newbedtool

    @_help('windowBed')
    @_file_or_bedtool()
    @_implicit('-a')
    def window(self,other, **kwargs):
        """
        Intersect with a window.

        Example usage::

            a = bedtool('in.bed')

            # Consider features up to 500 bp away as overlaps
            b = a.window(w=500)
        """
        tmp = self._tmp()
        cmds = ['windowBed',]
        kwargs['a'] = self.fn
        if (type(other) is str) or (type(other) is unicode):
            kwargs['b'] = other
        else:
            assert isinstance(other,bedtool), 'Either filename or another bedtool instance required'
            kwargs['b'] = other.fn
        cmds.extend(self.parse_kwargs(**kwargs))
        cmds.extend(['>',tmp])
        os.system(' '.join(cmds))
        return bedtool(tmp)

    @_help('shuffleBed')
    @_implicit('-i')
    def shuffle(self,genome=None,**kwargs):
        if genome is not None:
            genome_fn = self.get_genome(genome)
            kwargs['g'] = genome_fn
        kwargs['i'] = self.fn
        tmp = self._tmp()
        cmds = ['shuffleBed',]
        cmds.extend(self.parse_kwargs(**kwargs))
        cmds.extend(['>',tmp])
        os.system(' '.join(cmds))
        return bedtool(tmp)
    
    @_help('sortBed')
    @_implicit('-i')
    def sort(self,**kwargs):
        kwargs['i'] = self.fn
        cmds = ['sortBed']
        cmds.extend(self.parse_kwargs(**kwargs))
        tmp = self._tmp()
        cmds.extend(['>',tmp])
        os.system(' '.join(cmds))
        return bedtool(tmp)
    
    def features(self):
        """
        Returns an iterator of :class:`bedfeature` objects.
        """
        for line in self.__iterator():
            L = line.rstrip().split('\t')
            args = [None for i in range(12)]
            args[:len(L)] = L
            yield bedfeature(*args)

    def count(self):
        """
        Number of features in BED file. Does the same thing as len(self), which
        actually just calls this method.

        Only counts the actual features.  Ignores any track lines, browser
        lines, lines starting with a "#", or blank lines.

        Example usage::

            a = bedtool('in.bed')
            a.count()
        """
        c = 0
        f = open(self.fn)
        for i in f:
            if i.startswith('browser'):
                continue
            if i.startswith('track'):
                continue
            c += 1
        f.close()
        return c

    def print_sequence(self):
        """
        Print the sequence that was retrieved by the :meth:`bedtool.sequence`
        method.
        
        See usage example in :meth:`bedtool.sequence`.
        """
        if not hasattr(self,'seqfn'):
            raise ValueError, 'Use .sequence(fasta_fn) to get the sequence first'
        f = open(self.seqfn)
        s = f.read()
        f.close()
        return s

    def save_seqs(self,fn):
        """
        Save sequences of features in this bedtool object as a fasta file *fn*.

        In order to use this function, you need to have called
        the :meth:`bedtool.sequence()` method.
        
        A new bedtool object is returned which references the newly saved file.

        Example usage::

            a = bedtool('in.bed')

            # specify the filename of the genome in fasta format
            a.sequence('data/genomes/genome.fa')

            # use this method to save the seqs that correspond to the features
            # in "a"
            a.save_seqs('seqs.fa')
        """
        if not hasattr(self,'seqfn'):
            raise ValueError, 'Use .sequence(fasta_fn) to get the sequence first'
        fout = open(fn,'w')
        fout.write(open(self.seqfn).read())
        fout.close()
        return bedtool(fn)

    def get_genome(self, genome, fn=None):
        """
        Download chrom size info for *genome* from UCSC, removes the header
        line, and saves in a temp file.  Could be useful for end users,
        but mostly called internally by :meth:`bedtool.slop` and other methods that need
        the genome file.

        If *fn* is None, then saves as a temp file.

        Returns the filename.

        Example usage::

            a = bedtool('in.bed')
            fn = a.get_genome('dm3', 'dm3.genome')

        """
        tmp = self._tmp() 
        cmds = ['mysql',
                '--user=genome',
                '--host=genome-mysql.cse.ucsc.edu',
                '-A',
                '-e',
                '"select chrom, size from %s.chromInfo"' % genome,
                '|',
                'tail -n +2',
                '>',
                tmp]
        os.system(' '.join(cmds))
        self.genome_fn = tmp
        return tmp

    def newshuffle(self):
        """
        Quite fast implementation of shuffleBed; assumes 'dm3' as the genome
        and assumes shuffling within chroms.

        Example usage::

            a = bedtool('in.bed')

            # randomly shuffled version of "a"
            b = a.newshuffle()
        
        This is equivalent to the following command-line usage of ``shuffleBed``::

            shuffleBed -i in.bed -g dm3.genome -chrom -seed $RANDOM > /tmp/tmpfile

        """
        tmp = self._tmp()

        TMP = open(tmp,'w')
        for line in self.__iterator():
            chrom,start,stop = line.split()[:3]
            start = int(start)
            stop = int(stop)
            length = stop-start
            newstart = random.randint(self.genome[chrom][0], self.genome[chrom][1]-length)
            newstop = newstart + length
            TMP.write('%s\t%s\t%s\n' % (chrom,newstart,newstop))
        TMP.close()
        return bedtool(tmp)
 
    def randomstats(self, other, iterations, intersectkwargs={}):
        """
        Sends args to :meth:`bedtool.randomintersection` and compiles results
        into a dictionary with useful stats.  Requires scipy and numpy.

        Example usage::

            a = bedtool('in.bed')

            # Randomization results from 100 iterations, using the u=True kwarg (report
            # features in "a" only once for each intersection).
            results = a.randomstats('other.bed', iterations=100, intersectkwargs={'u':True})
        """
        try:
            from scipy import stats
            import numpy as np
        except ImportError:
            raise ImportError, "Need to install NumPy and SciPy for stats..."
        
        if (type(other) is str) or (type(other) is unicode):
            other = bedtool(other)
        else:
            assert isinstance(other,bedtool), 'Either filename or another bedtool instance required'

        counts = self.randomintersection(other,iterations,intersectkwargs)

        actual = len(self.intersect(other,**intersectkwargs))

        distribution = self.randomintersection(other, iterations=iterations, intersectkwargs=intersectkwargs)
        distribution = np.array(distribution)
        
        med_count = np.median(distribution)
        
        n = float(len(distribution))
        
        frac_above = sum(distribution >= actual)/n
        frac_below = sum(distribution <= actual)/n
        
        normalized = actual/med_count
        
        lower_95th = stats.scoreatpercentile(distribution,2.5)
        upper_95th = stats.scoreatpercentile(distribution,97.5)
        
        actual_percentile = stats.percentileofscore(distribution,actual)
        d = {
        'iterations':iterations,
        'actual': actual,
        'file_a':self.fn,
        'file_b':other.fn,
        self.fn: len(self),
        other.fn: len(other),
        'self':len(self),
        'other':len(other),
        'frac randomized above actual': frac_above,
        'frac randomized below actual': frac_below,
        'median randomized': med_count,
        'normalized': normalized,
        'lower_95th': lower_95th,
        'upper_95th': upper_95th,
        'percentile': actual_percentile,
        }
        return d

    def print_randomstats(self,other,iterations,intersectkwargs={}):
        """
        Nicely prints the reciprocal randomization of two files.
        """
        if (type(other) is str) or (type(other) is unicode):
            other = bedtool(other)
        d1 = self.randomstats(other, iterations, intersectkwargs)
        d2 = other.randomstats(self, iterations, intersectkwargs)
        
        print

        print 'Randomizing %s:' % self.fn
        print '\t%s features in %s' % (d1[self.fn],self.fn)
        print '\t%s features in %s' % (d1[other.fn],other.fn)
        print '\t%s actual intersections' % d1['actual']
        print '\t%.2f median randomized' % d1['median randomized']
        print '\t%.2f enrichment score' % d1['normalized']
        print '\t%.2f percentile' % d1['percentile']

        print 
        
        print 'Randomizing %s:' % other.fn
        print '\t%s features in %s' % (d2[other.fn],other.fn)
        print '\t%s features in %s' % (d2[self.fn],self.fn)
        print '\t%s actual intersection count' % d2['actual']
        print '\t%.2f median randomized' % d2['median randomized']
        print '\t%.2f enrichment score' % d2['normalized']
        print '\t%.2f percentile' % d2['percentile']

    def randomintersection(self, other, iterations, intersectkwargs={}):
        """
        Performs *iterations* shufflings of self, each time intersecting with
        *other*.  *intersectkwargs* are passed to self.intersect().  Returns a
        list of integers where each integer is the number of intersections of
        one shuffled file with *other*.
        
        Example usage::

            r = bedtool('in.bed').randomintersection('other.bed', 100, {'u': True})
        """
        counts = []
        for i in range(iterations):
            tmp = self.newshuffle()
            tmp2 = tmp.intersect(other,**intersectkwargs)
            counts.append(len(tmp2))
            os.remove(tmp.fn)
            os.remove(tmp2.fn)
            del(tmp)
            del(tmp2)
        return counts

    def deprecated_randomintersection(self,other,iterations,intersectkwargs={}):
        # list of intersection counts, one for each iteration
        if isinstance(other, basestring):
            other = bedtool(other)
        counts = []
        count1 = self.count()
        count2 = other.count()
        total = float(count1+count2)
        prob1 = count1/total
        for iteration in range(iterations):
            tmp1 = open(self._tmp(),'w') 
            tmp2 = open(self._tmp(),'w')
            for i in itertools.chain(self,other):
                if random.random() <= prob1:
                    tmp1.write(i)
                else:
                    tmp2.write(i)
            tmp1.close()
            tmp2.close()
            tmp_bed1 = bedtool(tmp1.name)
            tmp_bed2 = bedtool(tmp2.name)
            counts.append(len(tmp_bed1.intersect(tmp_bed2, **intersectkwargs)))
            os.remove(tmp1.name)
            os.remove(tmp2.name)

        return counts
                 

    @_file_or_bedtool()
    @_returns_bedtool()
    def cat(self,other, postmerge=True, **kwargs):
        """
        Concatenates two bedtools objects (or an object and a file) and does an
        optional post-merge of the features.

        Use *postmerge=False* if you want to keep features separate.

        TODO:

            currently truncates at BED3 format!

        kwargs are sent to :meth:`bedtool.merge`.

        Example usage::

            a = bedtool('in.bed')
            
            # concatenate and merge features together if they overlap and are
            # on the same strand
            b = a.cat('other.bed', s=True)
        """
        tmp = self._tmp()
        if (type(other) is str) or (type(other) is unicode):
            other = bedtool(other)
        else:
            assert isinstance(other,bedtool), 'Either filename or another bedtool instance required'
        TMP = open(tmp,'w')
        for line in self.__iterator():
            newline = '\t'.join(line.split()[:3])+'\n'
            TMP.write(newline)
        for line in other.__iterator():
            newline = '\t'.join(line.split()[:3])+'\n'
            TMP.write(newline)
        TMP.close()
        c = bedtool(tmp)
        if postmerge:
            d = c.merge(**kwargs)
            return d
        else:
            return c

            

    def tostring(self):
        '''
        Returns the BED file as a string.  You can also ``print`` the bedtool object
        to view its contents.

        Example usage::

            a = bedtool('in.bed')
            
            # this is one looong string which contains the entire file
            long_string = a.tostring()
        '''
        f = open(self.fn)
        s = f.read()
        f.close()
        return s

    @_returns_bedtool()
    def saveas(self,fn,trackline=None):
        """
        Save BED file as a new file, adding the optional *trackline* to the
        beginning.  
        
        Returns a new bedtool for the newly saved file.
        
        A newline is automatically added to the trackline if it does not
        already have one.

        Example usage::

            a = bedtool('in.bed')
            b = a.random_subset(5)
            b.saveas('random-5.bed',trackline='track name="random subset" color=128,128,255')
        """
        fout = open(fn,'w')
        if trackline is not None:
            fout.write(trackline.strip()+'\n')
        fout.write(self.tostring())
        fout.close()
        return bedtool(fn)

    @_file_or_bedtool()
    def intersection_report(self, other, basename=True, **kwargs):
        """
        Prints a report of the reciprocal intersections with another bed file
        or :class:`bedtool` object.

        If *basename* is True (default), only prints the basename of the file
        and not the whole path.

        a = bedtool('in.bed')
        a.intersection_report('other.bed')
        """
        if (type(other) is str) or (type(other) is unicode):
            other = bedtool(other)

        int1 = self.intersect(other, **kwargs).count()
        int2 = other.intersect(self.fn, **kwargs).count()
        
        count1 = self.count()
        count2 = other.count()
        
        self_fn = self.fn 
        other_fn = other.fn

        if basename:
            self_fn = os.path.basename(self_fn)
            other_fn = os.path.basename(other_fn)
         
        print '%s\n\t%s total\n\t%s (%.1f%%) of these intersect %s' % (self_fn, count1,  int1,  (float(int1)/count1)*100, other_fn)
        print '%s\n\t%s total\n\t%s (%.1f%%) of these intersect %s' % (other_fn, count2,  int2, (float(int2)/count2)*100, self_fn)

    @_returns_bedtool()
    def random_subset(self,n):
        '''
        Returns a new bedtools object containing a random subset of the
        features in this subset.  Currently does so by reading in all features;
        future updates should fix this to something more robust (e.g., newlines
        in a memory map)

        Example usage:: 

            a = bedtool('in.bed')
            
            # Choose 5 random features from 'in.bed'
            b = a.random_subset(5)
            
        '''
        features = list(self.__iterator())
        tmpfn = self._tmp()
        tmp = open(tmpfn,'w')
        for i in range(n):
            tmp.write(random.choice(features))
        tmp.close()
        return bedtool(tmpfn)
        

    def size_filter(self,min=0,max=1e15):
        """
        Returns a new bedtool object containing only those features that are 
        > *min* and < *max*.

        Example usage::

            a = bedtool('in.bed')

            # Only return features that are over 10 bp.
            b = a.size_filter(min=10)

        """
        tmpfn = self._tmp()
        tmp = open(tmpfn,'w')
        for feature in self.features():
            if min < len(feature) < max:
                tmp.write(feature.tostring())
        tmp.close()
        return bedtool(tmpfn)

    def sorted(self,col, reverse=None):
        '''Returns a new bedtool object, sorted by the column specified. col
        can be a list of columns.  BED columns that are ints (start, stop and
        value) will be sorted numerically; other columns will be
        alphabetical.
        
        reverse is a list of booleans, same length as col, specifying which 
        fields to reverse-sort.

        TODO: currently multiple columns aren't working!

        a = bedtool('in.fn')
        b = a.sorted(col=2) # sort by start position
        c = a.sorted(col=5,reverse=True) # reverse sort on the values
        '''
        
        if type(col) is not list:
            col = [col]

        if reverse is None:
            reverse = [False for i in col]
        elif type(reverse) is not list:
            reverse = [reverse]
        
        assert len(reverse) == len(col), 'reverse must be same length as col'
        
        if len(col) > 1:
            raise NotImplementedError,'multi-column sort not yet working correctly'

        d = {1:'1,1',
             2:'2n,2n',
             3:'3n,3n',
             4:'4,4',
             5:'5n,5n'}
        
        tmp = self._tmp()
        cmds = ['sort']
        for c,r in zip(col,reverse):
            if r:
                cmds.append('-k '+d[c]+'r')
            else:
                cmds.append('-k '+d[c])
        cmds.append(self.fn)
        cmds.extend( ['>',tmp] )
        os.system(' '.join(cmds))
        return bedtool(tmp)

    def sequence_coverage(self):
        """
        Returns the number of bases covered by this BED file.  Does a self.merge() first
        to remove potentially multiple-counting bases.

        Example usage::

            a = bedtool('in.bed')
            
            # total bp in genome covered by 'in.bed'
            total_bp = a.sequence_coverage()
        """
        b = self.merge()
        total_bp = 0
        for line in b.features():
            total_bp += len(feature)
        return total_bp

    def parse_kwargs(self,**kwargs):
        """Given a set of keyword arguments, turns them into a command
        line-ready list of strings.  E.g., the kwarg dict::

            kwargs = dict(c=True,f=0.5)

        will be returned as::

            ['-c','-f','0.5']

        If there are symbols (e.g., "|"), then the parameter is quoted."
        """
        illegal_chars = '!@#$%^&*(),-;:.<>?/|[]{} \'\\\"'
        cmds = []
        for key,value in kwargs.items():
            if value is True:
                cmds.append('-'+key)
                continue
            if (type(value) is tuple) or (type(value) is list):
                value = ','.join(map(str,value))
            if type(value) is str:
                for i in illegal_chars:
                    if i in value:
                        value = '"%s"' % value
                        break
            cmds.append('-'+key)
            cmds.append(str(value))
        return cmds

    @_returns_bedtool()
    def feature_centers(self,n,report_smaller=True):
        '''
        Returns a new bedtools object with just the centers of size n extracted
        from this object's features.

        If *report_smaller* is True, then report features that are smaller than
        *n*.  Otherwise, ignore them.

        Example usage::

            a = bedtool('in.bed')

            # 5bp on either side of the center of each feature
            b = a.feature_centers(100)
        '''
        tmpfn = self._tmp()
        tmp = open(tmpfn,'w')
        for line in self.__iterator():
            L = line.strip().split('\t')
            chrom,start,stop = L[:3]
            start = int(start)
            stop = int(stop)
            
            # if smaller than window size, decide whether to report it or not.
            if (stop-start) < n:
                if report_smaller:
                    tmp.write(line)
                    continue
                else:
                    continue

            left = floor(n/2.0)
            right = ceil(n/2.0)
            midpoint = start + (stop-start)/2
            newstart = str( int(midpoint - left))
            newstop = str( int(midpoint + right))
            L[1] = newstart
            L[2] = newstop
            tmp.write('\t'.join(L)+'\n')
        tmp.close()
        return bedtool(tmpfn)

    @_returns_bedtool()
    def rename_features(self,new_name):
        """
        Forces a rename of all features.  Useful for if you have a BED file of
        exons and you want all of them to have the name "exon".
        """
        tmpfn = self._tmp()
        tmp = open(tmpfn,'w')
        for line in self.__iterator():
            L = line.strip().split('\t')
            chrom,start,stop = L[:3]
            if len(L) > 3:
                L[3] = new_name
            else:
                L.append(new_name)
            tmp.write('\t'.join(L)+'\n')
        tmp.close()
        return bedtool(tmpfn)
    
    @_returns_bedtool()
    def with_attrs(self, **kwargs):
        """
        Given arbitrary keyword arguments, turns the keys and values into
        attributes.

        Example usage::

            # add a "label" attribute to each bedtool
            a = bedtool('a.bed').with_attrs(label='transcription factor 1')
            b = bedtool('b.bed').with_attrs(label='transcription factor 2')
            for i in [a,b]:
                print i.count(), 'features for', i.label
        """
        for key,value in kwargs.items():
            setattr(self,key,value)
        return self

    def counts(self):
        """
        After running :meth:`bedtool.intersect` with the kwarg *c=True*, use
        this method to return a list of the count of features in "b" that
        intersected each feature in "a".

        Example usage::

            a = bedtool('in.bed')
            b = a.intersect('other.bed', c=True)
            counts = b.counts()

            # assuming you have matplotlib installed, plot a histogram
            
            import pylab
            pylab.hist(counts)
            pylab.show()
        """
        if not self._hascounts:
            raise ValueError, 'Need intersection counts; run intersection(fn, c=True) for this or manually set self._hascounts=True.'
        counts = []
        for line in self.__iterator():
            L = line.split()
            chrom,start,stop = L[:3]
            count = int(L[-1])
            counts.append(count)
        return counts

    def normalized_counts(self):
        """
        After running :meth:`bedtool.intersect` with the kwarg *c=True*, use
        this method to return a list of the density of features in "b" that
        intersected each feature in "a".

        This takes the counts in each feature and divides by the bp in that
        feature.

        Example usage::

            a = bedtool('in.bed')
            b = a.intersect('other.bed', c=True)
            counts = b.normalized_counts()

            # assuming you have matplotlib installed, plot a histogram
            
            import pylab
            pylab.hist(counts)
            pylab.show()
        """
        if not self._hascounts:
            raise ValueError, 'Need intersection counts; run intersection(fn, c=True) for this or manually set self._hascounts=True.'
        normalized_counts = []
        for line in self.__iterator():
            L = line.split()
            chrom,start,stop = L[:3]
            count = float(L[-1])
            normalized_count = count/(int(stop)-int(start))*1000
            normalized_counts.append(normalized_count)
        return normalized_counts

    def lengths(self):
        """
        Returns a list of feature lengths.

        Example usage::
            
            a = bedtool('in.bed')

            lengths = a.lengths()

            # if you have pylab installed, plot a histogram
            import pylab
            pylab.hist(lengths)
            pylab.show()
        """
        feature_lengths = []
        for line in self.__iterator():
            chrom,start,stop = line.split()[:3]
            start = int(start)
            stop = int(stop)
            length = stop-start
            feature_lengths.append(length)
        return feature_lengths
    
