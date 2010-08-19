"""
pybedtools wraps Aaron Quinlan's BEDtools programs in an easy-to-use class for
manipulating BED files from within Python.

The bedtool class includes some functionality not included in BEDtools (size
filters, random selection of features, automatic chromosome-file retrieval,
etc)

Ryan Dale
"""
import sys
import os
import tempfile
import subprocess
import random

# Module-level list of all tempfiles created.  These will be deleted when
# cleanup() is called.
TEMPFILES = []


def set_tempdir(tempdir):
    """
    Sets the directory for temp files.  Useful for clusters that use a /scratch
    partition rather than a /tmp dir.  Convenience function to simply set
    tempfile.tempdir.
    """
    tempfile.tempdir = tempdir

def cleanup():
    """Deletes all temporary files in *TEMPFILES*"""
    for fn in TEMPFILES:
        if os.path.exists(fn):
            os.remove(fn)

def help(command):
    '''Decorator that adds help from each of the BEDtools programs to the
    docstring of the method that calls the program'''
    p = subprocess.Popen([command], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    help = p.communicate()[1]
    help = help.replace('_','**')
    # insert tabs into the help
    help = help.split('\n')
    help = ['\t'+i for i in help]
    help = '\n'.join(help)
    def decorator(func):
        if func.__doc__ is None:
            func.__doc__ = ''
        orig = func.__doc__
        func.__doc__ = '*pybedtools help:*\n'
        func.__doc__ += orig
        func.__doc__ += '\n\n*Original BEDtools program help:*\n'
        func.__doc__ += help
        return func
    return decorator

class bedfeature(object):
    def __init__(self, chr,start,stop,
                 name=None,value=None,strand=None,
                 thickStart=None,thickStop=None,itemRGB=None,
                 blockCount=None,blockSizes=None,blockStarts=None):
        self.chr=chr
        self.start=int(start)
        self.stop=int(stop)
        self.name=name
        self.strand=strand
        try:
            self.value=float(value)
        except TypeError:
            self.value=value
        try:
            self.thickStart=int(thickStart)
        except TypeError:
            self.thickStart=thickStart
        try:
            self.thickStop=int(thickStop)
        except TypeError:
            self.thickStop=thickStop
        try:
            self.blockCount=int(blockCount)
        except TypeError:
            self.blockCount=blockCount
        
        self.itemRGB=itemRGB
        self.blockSizes=blockSizes
        self.blockStarts=blockStarts

    def __repr__(self):
        return 'bed feature: %s:%s-%s' % (self.chr,self.start,self.stop)
    
    def __len__(self):
        return self.stop - self.start

    def tostring(self):
        """Prints the bed record suitable for writing to file, newline included.
        
        In the interest of speed, does not do error-checking.
        """
        items = [self.chr, 
                 self.start, 
                 self.stop, 
                 self.name, 
                 self.value, 
                 self.strand, 
                 self.thickStart,
                 self.thickStop, 
                 self.itemRGB,
                 self.blockCount,
                 self.blockSizes, 
                 self.blockStarts]
        printables = []
        for item in items:
            if item is None:
                printables.append('')
            else:
                printables.append(str(item))
        
        return '\t'.join(printables).rstrip()+'\n'


class bedtool(object):
    """
    Wrapper around ``BEDtools`` suite of programs; also contains many useful
    methods for more detailed work with BED files.

    Example usage::

        a = bedtool('in.bed')
    """
    def __init__(self,fn):
        """
        *fn* is a BED format file.
        """
        if not os.path.exists(fn):
            raise ValueError, 'File "%s" does not exist' % fn
        self.fn = fn
        self._hascounts = False

    def _tmp(self):
        '''
        Makes a tempfile and registers it for eventual deletion when the
        instance is deleted.  Adds a "pybedtools." prefix and ".tmp" extension
        for easy deletion if you forget to call pybedtools.cleanup().
        '''

        tmpfn = tempfile.NamedTemporaryFile(prefix='pybedtools.',suffix='.tmp',delete=False)
        tmpfn = tmpfn.name
        TEMPFILES.append(tmpfn)
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

    def __repr__(self):
        return '<bedtool (%s)>'%self.fn

    def __str__(self):
        f = open(self.fn)
        s = f.read()
        f.close()
        return s

    def __len__(self):
        return self.count()

    @property
    def lines(self):
        return open(self.fn)
   
    def head(self,n=10):
        for i,line in enumerate(self.lines):
            if i == (n):
                break
            print line,

    @help('intersectBed')
    def intersect(self, other, **kwargs):
        """
        Intersect with another BED file. If you want to use BAM, specify
        *abam='filename.bam'*.  Returns a new bedtool object.
       
        Example usage::

            # create bedtool object
            a = bedtool('in.bed')

            # get overlaps with other.bed
            overlaps = a.intersect('other.bed')

            # use v=True to get the inverse, or those unique to in.bed
            unique_to_a = a.intersect('other.bed', v=True)

            # create a new bedtool object and intersect it with a
            # to get the feature unique to this bed file
            unique_to_other = bedtool('other.bed').intersect(a, v=True)

        """
        if type(other) is str:
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

    @help('fastaFromBed')
    def sequence(self, **kwargs):
        '''
        Wraps ``fastaFromBed``.  *fi* is passed in by the user; *bed* is
        automatically passed in as the bedfile of this object; *fo* by default
        is a temp file.  Use save_seqs() to save as a file.
        
        Example usage::

            a = bedtool('in.bed')
            a.sequence('genome.fa')
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

    @help('subtractBed')
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
        if type(other) is str:
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

    @help('slopBed')
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

    @help('mergeBed')
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

    @help('closestBed')
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
        if type(other) is str:
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

    @help('windowBed')
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
        if type(other) is str:
            kwargs['b'] = other
        else:
            assert isinstance(other,bedtool), 'Either filename or another bedtool instance required'
            kwargs['b'] = other.fn
        cmds.extend(self.parse_kwargs(**kwargs))
        cmds.extend(['>',tmp])
        os.system(' '.join(cmds))
        return bedtool(tmp)

    @help('groupBy')
    def groupBy(self,**kwargs):
        tmp = self._tmp()
        cmds = ['groupBy']
        kwargs['i'] = self.fn
        cmds.extend(self.parse_kwargs(**kwargs))
        cmds.extend(['>',tmp])
        os.system(' '.join(cmds))
        return bedtool(tmp)

    @help('shuffleBed')
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
    
    @help('sortBed')
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

    def print_sequence(self,fn):
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

    def get_genome(self,genome):
        """
        Download chrom size info for *genome* from UCSC, removes the header
        line, and saves in a temp file.  Could be useful for end users,
        but mostly called internally by :meth:`bedtool.slop` and other methods that need
        the genome file.

        Example usage::

            a = bedtool('in.bed')
            fn = a.get_genome('dm3')

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
        dm3 = """
        chr2L     23011544
        chr2LHet  368872
        chr2R     21146708
        chr2RHet  3288761
        chr3L     24543557
        chr3LHet  2555491
        chr3R     27905053
        chr3RHet  2517507
        chr4      1351857
        chrU      10049037
        chrUextra 29004656
        chrX      22422827
        chrXHet   204112
        chrYHet   347038
        chrM      19517
        """
        tmp = self._tmp()
        chromlens = {} 
        for i in dm3.split('\n'):
            i = i.strip()
            if len(i) < 1:
                continue
            chrom,length = i.split()
            chromlens[chrom] = int(length)

        TMP = open(tmp,'w')
        for line in self.__iterator():
            chrom,start,stop = line.split()[:3]
            start = int(start)
            stop = int(stop)
            length = stop-start
            newstart = random.randint(1, chromlens[chrom]-length)
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
        from scipy import stats
        import numpy as np
        
        if type(other) is str:
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
        'actual': actual,
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
        if type(other) is str:
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
        if type(other) is str:
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

    def intersection_report(self, other, basename=True, **kwargs):
        """
        Prints a report of the reciprocal intersections with another bed file
        or :class:`bedtool` object.

        If *basename* is True (default), only prints the basename of the file
        and not the whole path.

        a = bedtool('in.bed')
        a.intersection_report('other.bed')
        """
        if type(other) is str:
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
                tmp.write(line)
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
        """
        cmds = []
        for key,value in kwargs.items():
            if value is True:
                cmds.append('-'+key)
                continue
            if (type(value) is tuple) or (type(value) is list):
                value = ','.join(map(str,value))
            cmds.append('-'+key)
            cmds.append(str(value))
        return cmds

    def feature_centers(self,n):
        '''
        Returns a new bedtools object with just the centers of size n extracted from
        this object's features.

        Example usage::

            a = bedtool('in.bed')

            # 100bp on either side of the center of each feature
            b = a.feature_centers(100)
        '''
        tmpfn = self._tmp()
        tmp = open(tmpfn,'w')
        for line in self.__iterator():
            L = line.strip().split('\t')
            chrom,start,stop = L[:3]
            start = int(start)
            stop = int(stop)
            halfwidth = (stop-start)/2
            if halfwidth < 2*n:
                continue
            midpoint = start+halfwidth
            newstart = str(midpoint - n)
            newstop = str(midpoint + n)
            L[1] = newstart
            L[2] = newstop
            tmp.write('\t'.join(L)+'\n')
        tmp.close()
        return bedtool(tmpfn)

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
            chrom,start,stop,count = line.split()
            count = int(count)
            counts.append(count)
        return counts

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

