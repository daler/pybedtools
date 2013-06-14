import os
import pybedtools
import itertools
from collections import defaultdict


class MultiClassifier(object):
    def __init__(self, bed, annotations, genome, sample_name='sample', names=None,
                 prefix='.'):
        """
        Classifies files using bedtools multiinter.

        An example use case is to classify reads as intronic or exonic.

        `bed` must be a BED/VCF/GFF/GTF BedTool object or filename -- not
        a BAM.  If you want to use reads from a BAM file, first convert to
        bed::

            bed = pybedtools.BedTool('reads.bam').bam_to_bed(split=True).fn

        `annotations` will be used to classify the features in `bed`.

            If `annotations` is a list, then use those annotations directly.

            If `annotations` is a BedTool or string, then split the annotations
            into multiple files based on the featuretype (assumes GFF/GTF
            format).  Each of these new filenames will get the `prefix`
            supplied.

        `sample_name` is the name of the sample -- used internally, but it
        needs to be a unique name that is not the name of a featuretype in the
        annotations file (i.e., sample_name="exon" would not be a good choice).

        `names` is an optional list of names that correspond to `annotations`
        if it is a list, OR a function that maps items in `annotations` to
        a featuretype name -- e.g.::

            names=lambda x: x.replace('.gff', '')

        `genome` is the string assembly name, or dictionary of {chrom: (start,
        stop)} coordinates for each chromosome.  It is used to determine the
        "unannotated" space of the genome.

        After running the `classify()` method, the results are available as
        self.results, and these are parsed into two dictionaries,
        self.class_genome_bp, which holds the total number of annotated genomic
        bp for each class, and self.class_sample_bp, which holds the total
        number of bp in the sample that overlapped each class.

        The table() method offers a means to filter/merge these
        fully-classified dictionaries such that they only contain featuretypes
        of interest (by using the `include` kwarg).

        The print_table() method prints a nice, optionally-filtered/merged
        table of the results, sorted by total sample bp.

        Example usage::

            >>> bam = pybedtools.example_bedtool('x.bam')
            >>> bed = bam.bam_to_bed(split=True).sort()
            >>> anno = pybedtools.example_filename('dm3-chr2L-5M.gff.gz')
            >>> names, fns = MultiClassifier.split_annotations(anno)

            >>> # Example method of making a name converter.
            >>> # (note that the split `fns` have the form 'split_exon.gff')
            >>> def name_converter(x):
            ...     return x.split('_')[-1]

            >>> c = MultiClassifier(bed, fns, names=name_converter, genome='dm3')
            >>> c.classify()
            >>> inc = ['exon', 'intron']
            >>> table = c.print_table(include=inc)
            >>> print table #doctest: +NORMALIZE_WHITESPACE
            class	sample	genome
            exon	337580	618147
            unannotated	16424	18105002
            exon, intron	12882	118826
            intron	9949	375417

            >>> # Clean up the split files.
            >>> for fn in fns:
            ...     os.unlink(fn)

        """
        self.bed = pybedtools.BedTool(bed)
        self.sample_name = sample_name
        self.genome = genome

        # If a list of annotation files was provided, then use them directly,
        # ensuring that they are BedTools objects
        if isinstance(annotations, (list, tuple)):
            annotations = [pybedtools.BedTool(i).fn for i in annotations]

            # ... name-munge as necessary
            if names is None:
                names = annotations
            if hasattr(names, '__call__'):
                names = [names(i) for i in annotations]

        # Otherwise, split the single annotation file into a form suitible for
        # use by multiintersect
        else:
            names, files = self.split_annotations(annotations)
            annotations = list(files)
            names = list(names)

        # multiintersect still needs the input file (e.g., peaks)
        self.annotations = [self.bed.fn] + annotations
        self.names = [self.sample_name] + names

        self.class_sample_bp = defaultdict(int)
        self.class_genome_bp = defaultdict(int)

    @classmethod
    def split_annotations(self, annotations, prefix='split_', name_func=None):
        """
        Given an annotation file in GFF format, split into different files --
        each file containing only one type of feature.

        `prefix` will be added to each featuretype to construct the filenames.

        `name_func`, by default, will use the feature type field of a GTF/GFF
        file, or the feature.name attribute of another format.  Supply a custom
        function that accepts a pybedtools.Interval instance for more control.
        """
        # dict of open files, one for each featuretype found
        files = {}
        annotations = pybedtools.BedTool(annotations)
        ft = annotations.file_type
        if name_func is None:
            if ft == 'gff':
                name_func = lambda x: x[2]
            else:
                name_func = lambda x: x.name
        for feature in annotations:
            featuretype = name_func(feature)
            if featuretype not in files:
                filename = '%s%s' % (prefix, featuretype)
                files[featuretype] = open(filename, 'w')
            files[featuretype].write(str(feature))
        for f in files.values():
            f.close()
        return zip(*[(featuretype, f.name)
                   for (featuretype, f) in files.items()])

    def classify(self, **kwargs):
        """
        Classify the features in self.bed, populating several dictionaries in
        self.  `kwargs` are passed on to BedTool.multi_intersect.  The
        "empty=True" kwarg to multi_intersect is always provided to make sure
        the classification works correctly.
        """
        self.results = self.bed.multi_intersect(
            i=[self.bed.fn] + self.annotations,
            names=[self.sample_name] + self.names,
            genome=self.genome,
            empty=True)

        sample = set([self.sample_name])

        for i in self.results:
            # Even if there were multiple annotations, only report *that* there
            # was a hit, e.g., "exon,exon,exon,5'UTR" -> (exon, 5'UTR)
            full_class = frozenset(i[4].split(','))

            # Including sample name in class name would be redundant, so remove
            # it
            class_name = full_class.difference(sample)

            # Only report if sample was in the class
            if self.sample_name in full_class:
                self.class_sample_bp[class_name] += len(i)

            # But always report the presence of the class, regardless of if
            # there was a hit in the sample or not.
            self.class_genome_bp[class_name] += len(i)

        # Genomic unannotated has the key ["none"]; sample unannotated as the
        # key [] (because it was set-differenced out)
        self.class_genome_bp[frozenset(['unannotated'])] = \
            self.class_genome_bp.pop(frozenset(['none']), 0) \
            + self.class_genome_bp.pop(frozenset([]), 0)


        self.class_sample_bp[frozenset(['unannotated'])] = \
                             self.class_sample_bp.pop(frozenset([]), 0)

    def table(self, include=None):
        """
        If `include` is not None, then return versions of self.class_genome_bp
        and self.class_sample_bp that only look at the featuretypes in
        `include`.

        Otherwise, simply return these dictionaries unchanged (and including
        all available featuretypes)
        """
        if not include:
            return self.class_sample_bp, self.class_genome_bp

        d = defaultdict(int)
        s = defaultdict(int)
        include = set(include)
        for key in self.class_genome_bp:
            # orig data
            seq_bp = self.class_genome_bp[key]
            bed_bp = self.class_sample_bp[key]

            # create a merged class name by keeping only featuretypes in
            # `include`
            merged_class = frozenset(set(key).intersection(include))

            # update the new dictionaries
            d[merged_class] += bed_bp
            s[merged_class] += seq_bp

        return d, s

    def hierarchical_table(self, order, include=None, include_unannotated=True):
        """
        Returns a hierarchically-ordered table, using the specified `order`.

        For example the order::

                ['TSS', 'five_prime_UTR', 'CDS', 'intron', 'three_prime_UTR', 'TTS']

        would weight the classes from the 5'-most end of the gene.

        This summarizes the classes based on the highest-priority featuretype
        in the hierarchy.

        For example, using the above hierarchy, the following summary-class
        assignments will be made::

            (TSS, five_prime_UTR)          -> TSS
            (intron, CDS, TSS)             -> TSS
            (intron, CDS, three_prime_UTR) -> CDS

        The table has the following format, where the "classes" list is ordered by sample bp.

            {
                'TSS': [
                        (<class name 1>, <sample bp>, <genomic bp>),
                        (<class name 2>, <sample bp>, <genomic bp>),
                        ...
                        ],

                'five_prime_UTR': [
                        (<class name 1>, <sample bp>, <genomic bp>),
                        (<class name 2>, <sample bp>, <genomic bp>),
                        ...
                        ],
                ...

            }
        """
        sample, genomic = self.table(include=include)
        sample = sample.copy()
        genomic = genomic.copy()
        keys = list(set(sample.keys() + genomic.keys()))
        table = {}
        for h in order:
            classes = []
            for k in keys:
                if h in k:
                    sample_bp = sample.pop(k, 0)
                    genomic_bp = genomic.pop(k, 0)
                    classes.append(
                        (k, sample_bp, genomic_bp)
                    )
            table[h] = sorted(classes, key=lambda x: x[1], reverse=True)
        #if include_unannotated:
        #    table['unannotated'] = [(frozenset(['unannotated']), sample.pop(frozenset([]), 0), genomic.pop(frozenset([]), 0))]
        return table


    def print_table(self, include=None):
        """
        Returns a string containing a tab-delimited table, including header,
        with classes sorted by total bp in each class.
        """
        d, s = self.table(include)
        out = []
        out.append('class\tsample\tgenome')
        for cls in sorted(d.keys(), key=lambda x: d[x], reverse=True):
            if len(cls) == 0:
                label = 'unannotated'
            else:
                label = ', '.join(sorted(cls))

            out.append('%s\t%s\t%s' % (label, d[cls], s[cls]))
        return '\n'.join(out)
