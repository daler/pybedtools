import os
import pybedtools
import itertools
from collections import defaultdict


class BasePairClassifier(object):
    def __init__(self, bed, annotations, genome, sample_name='sample',
                 names=None, prefix='split_'):
        """
        Classifies files using bedtools multiinter.

        The results from this class split up the bp in `bed` in to classes as
        annotated in `annotations`.  Note that this is different from counting
        the number of features in `bed` that fall within `annotations`; for
        this latter functionality, see the FeatureClassifier class in this same
        module.

        `bed` must be a BED/VCF/GFF/GTF BedTool object or filename -- not
        a BAM.  This is because `bedtools multiinter`, the program called here,
        does not accept BAM as input. If you want to use reads from a BAM file,
        first convert to bed, i.e.,::

            bed = pybedtools.BedTool('reads.bam').bam_to_bed(split=True).fn

        `annotations` will be used to classify the features in `bed`.

            If `annotations` is a list, then use those annotations directly.

            If `annotations` is a BedTool or string, then split the annotations
            into multiple files based on the featuretype (assumes GFF/GTF
            format).  Each of these new filenames will get the `prefix`
            supplied.

        `sample_name` is the name of the sample. It is used internally, but it
        needs to be a unique name that is not the name of a featuretype in the
        annotations file (i.e., sample_name="exon" would not be a good choice).

        `names` is an optional list of names that correspond to `annotations`
        if it is a list, OR a function that maps items in `annotations` to
        a featuretype name -- e.g.::

            names=lambda x: x.replace('.gff', '')

        `genome` is the string assembly name, or dictionary of {chrom: (start,
        stop)} coordinates for each chromosome.  It is used to determine the
        "unannotated" space of the genome.

        After running the `classify()` method, the results BedTool is available
        as self.results, and these are parsed into two dictionaries,
        1) self.class_genome_bp, which holds the total number of annotated
        genomic bp for each class, and 2) self.class_sample_bp, which holds the
        total number of bp in the sample that overlapped each class.

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
            >>> def f(x):
            ...     return x.split('_')[-1]

            >>> c = MultiClassifier(bed, fns, names=f, genome='dm3')
            >>> c.classify()
            >>> inc = ['exon', 'intron']
            >>> table = c.print_table(include=inc)
            >>> print table #doctest: +NORMALIZE_WHITESPACE
            class	sample	genome
            exon	361971	1188410
            unannotated	19103	19841604
            exon, intron	14250	177913
            intron	10657	1803617

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
        for f in list(files.values()):
            f.close()
        return list(zip(*[(featuretype, f.name)
                   for (featuretype, f) in list(files.items())]))

    def classify(self, **kwargs):
        """
        Classify the features in self.bed, populating several dictionaries in
        self.  `kwargs` are passed on to BedTool.multi_intersect.  The
        "empty=True" kwarg to multi_intersect is always provided to make sure
        the classification works correctly.
        """
        self.results = self.bed.multi_intersect(
            i=self.annotations,
            names=self.names,
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

    def hierarchical_table(self, order, include=None,
                           include_unannotated=True):
        """
        Returns a hierarchically-ordered table, using the specified `order`.

        For example the order::

                ['TSS', 'five_prime_UTR', 'CDS', 'intron', 'three_prime_UTR',
                'TTS']

        would weight the classes from the 5'-most end of the gene.

        This summarizes the classes based on the highest-priority featuretype
        in the hierarchy.

        For example, using the above hierarchy, the following summary-class
        assignments will be made::

            (TSS, five_prime_UTR)          -> TSS
            (intron, CDS, TSS)             -> TSS
            (intron, CDS, three_prime_UTR) -> CDS

        The table has the following format, where the "classes" list is ordered
        by sample bp.

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
        keys = list(set(list(sample.keys()) + list(genomic.keys())))
        table = {}
        for h in order:
            classes = []
            for k in keys:
                if h in k:
                    try:
                        sample_bp = sample.pop(k, 0)
                        genomic_bp = genomic.pop(k, 0)
                        classes.append(
                            (k, sample_bp, genomic_bp)
                        )
                    except KeyError:
                        # already popped
                        continue
            table[h] = sorted(classes, key=lambda x: x[1], reverse=True)
        if include_unannotated:
            table['unannotated'] = [(frozenset(['unannotated']),
                                     sum(sample.values()),
                                     sum(genomic.values()))]
        return table

    def print_table(self, include=None):
        """
        Returns a string containing a tab-delimited table, including header,
        with classes sorted by total bp in each class.
        """
        d, s = self.table(include)
        out = []
        out.append('class\tsample\tgenome')
        for cls in sorted(list(d.keys()), key=lambda x: d[x], reverse=True):
            if len(cls) == 0:
                label = 'unannotated'
            else:
                label = ', '.join(sorted(cls))

            out.append('%s\t%s\t%s' % (label, d[cls], s[cls]))
        return '\n'.join(out)

MultiClassifier = BasePairClassifier


class Classifier(object):
    """
    Classify intervals in one file by featuretypes in another.
    """
    def __init__(self, bed, annotations):
        """
        Classify features in `bed` -- typically a BED or SAM/BAM but can be any
        format supported by BedTools -- into classes based on the featuretypes
        in the GFF/GTF file, `annotations`.

        For example, you can classify ChIP-seq peaks in a BED file by intron,
        exon, or whatever is annotated in the GFF file.  If you want to
        consider promoter regions, you'll have to add these features yourself.

        The `class_counts` dictionary has its keys as sets of featuretypes
        (each one can be considered a "class" of features) and the value is the
        number of features in that class.  The special empty set class contains
        features that did not fall in an annotated region.

        You can access the individual features in the `class_features`
        dictionary, which contains the same keys but instead of counts, it
        contains the features themselves.  This is nice for saving the features
        in a separate BED file, e.g.,

        Furthermore, you can look up the class of any feature in the original
        BED file using the `feature_classes` dictionary::


        Example usage::

            >>> bed = pybedtools.example_filename('gdc.bed')
            >>> gff = pybedtools.example_filename('gdc.gff')
            >>> c = pybedtools.contrib.classifier.Classifier(bed, gff)
            >>> c.classify(include=['intron', 'exon'])
            >>> results = c.class_counts
            >>> results == {
            ... frozenset([]): 1,
            ... frozenset(['exon']): 3,
            ... frozenset(['intron']): 3,
            ... frozenset(['intron', 'exon']): 1}
            True
            >>> key = frozenset(['intron'])
            >>> features = c.class_features[key]
            >>> pybedtools.BedTool(iter(features)).saveas() #doctest: +ELLIPSIS
            <BedTool(...)>
            >>> feature = pybedtools.BedTool(bed)[2]
            >>> c.feature_classes[feature] == set(['intron', '.'])
            True

        """
        self.bed = pybedtools.BedTool(bed)
        self.annotations = pybedtools.BedTool(annotations)
        if self.annotations.file_type != 'gff':
            raise ValueError('Annotations file must be a GFF or GTF file; '
                             '%s appears to be a %s file' % (
                                 annotations,
                                 self.annotations.file_type))

    def available_featuretypes(self):
        """
        List the featuretypes available in the annotations file.
        """
        featuretypes = set()
        for i in self.annotations:
            featuretypes.update([i[2]])
        return list(featuretypes)

    def classify(self, include=None, exclude=None, stranded=False):
        """
        Perform classification, populating dictionaries in `self`.

        Intersect the BED file with the annotations file and return
        a dictionary where keys are BED features and values are the set of
        featuretypes that BED feature was found in.


        `include` is an optional list of featuretypes to restrict the
        classification to

        `exclude` is an optional list of featuretypes to exclude from
        classification (all other featuretypes will be used).

        To see what's available, use available_featuretypes().

        When run, this method creates the following dictionaries as attributes
        of this object:

         :feature_classes:
            keys are Intervals from `bed`; values are sets of featuretypes from
            `annotations`

         :class_features:
            keys are frozensets of featuretypes from `annotations`; values are
            lists of Intervals from `bed`;

         :class_counts:
            keys are frozensets of featuretypes from annotations`; values are
            number of features -- so the length of values in the class_features
            dictionary.

        """
        if include and exclude:
            raise ValueError('Can only specify one of `include` or `exclude`')
        if exclude:
            exclude = set(exclude)
        if include:
            include = set(include)

        # Figure out the index of the featuretype field in the output
        bed_fields = self.bed.field_count()
        featuretype_idx = bed_fields + 2

        self.feature_classes = defaultdict(set)

        x = self.bed.intersect(self.annotations,
                               wao=True,
                               s=stranded,
                               stream=True)
        for feature in x:
            featuretype = feature[featuretype_idx]

            # If we're not supposed to consider this featuretype, then reset to
            # the standard GFF empty field string of "."
            if (include and featuretype not in include) \
                    or (exclude and featuretype in exclude):
                featuretype = '.'

            # the original query is in the first `bed_fields` items.  Construct
            # a new feature out of this and use it as the key.
            key = pybedtools.create_interval_from_list(
                feature.fields[:bed_fields])
            self.feature_classes[key].update([featuretype])

        self.class_features = defaultdict(list)
        self.class_counts = defaultdict(int)

        for feature, featuretypes in list(self.feature_classes.items()):
            # get rid of "unannotated"
            ft = featuretypes.difference(['.'])
            key = frozenset(ft)
            self.class_features[key].append(feature)
            self.class_counts[key] += 1

        # convert defaultdicts to regular dicts
        self.class_features = dict(self.class_features)
        self.class_counts = dict(self.class_counts)
        self.feature_classes = dict(self.feature_classes)

    def features_to_file(self, prefix="", suffix=""):
        """
        Writes a set of files, one for each class.

        The filenames will be constructed based on the class names, track lines
        will be added to indicate classes, and `prefix` and `suffix` will be
        added to the filenames.
        """
        def make_filename(klass):
            return prefix + '_'.join(sorted(list(klass))) + suffix

        def make_trackline(klass):
            return 'track name="%s"' % (' '.join(sorted(list(klass))))

        for klass, features in self.class_features.items():
            pybedtools.BedTool(features)\
                .saveas(make_filename(klass), make_trackline(klass))

    def hierarchical_table(self, order, include=None,
                           include_unannotated=True):
        """
        Returns a hierarchically-ordered table, using the specified `order`.

        For example the order::

            ['TSS', 'five_prime_UTR', 'CDS', 'intron', 'three_prime_UTR',
            'TTS']

        would weight the classes from the 5'-most end of the gene.

        This summarizes the classes based on the highest-priority featuretype
        in the hierarchy.

        For example, using the above hierarchy, the following summary-class
        assignments will be made::

            (TSS, five_prime_UTR)          -> TSS
            (intron, CDS, TSS)             -> TSS
            (intron, CDS, three_prime_UTR) -> CDS

        The table has the following format, where the "classes" list is ordered
        by sample bp.

            {
                'TSS': [
                        (<class name 1>, count),
                        (<class name 2>, count),
                        ...
                        ],

                'five_prime_UTR': [
                        (<class name 1>, count),
                        (<class name 2>, count),
                        ...
                        ],
                ...

            }
        """
        counts = self.table(include=include)
        counts = counts.copy()
        keys = list(counts.keys())
        table = {}
        for h in order:
            classes = []
            for k in keys:
                if h in k:
                    try:
                        count = counts.pop(k)
                        classes.append((k, count))
                    except KeyError:
                        # i.e., already popped off
                        continue
            table[h] = sorted(classes, key=lambda x: x[1], reverse=True)
        if include_unannotated:
            table['unannotated'] = [(frozenset(['unannotated']), sum(counts.values()))]
        return table

    def table(self, include=None):
        """
        If `include` is not None, then return a copy of self.class_counts that
        only have at the featuretypes in `include`.

        Otherwise, return a simple copy.
        """
        if not include:
            return self.class_counts.copy()

        d = defaultdict(int)
        s = defaultdict(int)
        include = set(include)
        for key in self.class_counts:
            # orig data
            c = self.class_counts[key]

            # create a merged class name by keeping only featuretypes in
            # `include`.  If none of the members were in `include`, then the
            # counts are appended to entry keyed by set().
            merged_class = frozenset(set(key).intersection(include))

            # update the new dictionaries
            d[merged_class] += c

        return d


    def pie(self, hierarchical_table, ax=None, order=None, **kwargs):
        """
        `hierarchical_table` is the result of calling self.hierarchical_table()`.

        `ax` is an optional Axes object to plot onto, otherwise a new figure and axes will be created.

        `order`, if None, will sort the items from largest to smallest.
        Otherwise, `order` is the order of keys in `hierarchical_table` in
        which they should be plotted, counter-clockwise from the positive
        x axis.

        Additional args are passed to ax.pie.
        """
        from matplotlib import pyplot as plt

        # restructure the table so that it's a dict of class sums (rather than
        # having every class)
        d = {}
        for k, v in list(hierarchical_table.items()):
            d[k] = sum(i[1] for i in v)

        total = float(sum(d.values()))
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.axis('equal')

        if order:
            items = [(k, d[k]) for k in order]
        else:
            items = sorted(list(d.items()), key=lambda x: x[1])

        newlabels, labels, counts = [], [], []
        for label, count in items:
            labels.append(label)
            counts.append(count)
            frac = count / total * 100.
            newlabels.append('{label}: {count} ({frac:.1f}%)'.format(**locals()))
        ax.pie(
            x=counts, labels=newlabels, labeldistance=1.2, **kwargs)
        return ax
