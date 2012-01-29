import pybedtools
import itertools
from collections import defaultdict


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
            >>> c = pybedtools.contrib.Classifier(bed, gff)
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
            >>> pybedtools.BedTool(iter(features)).saveas()  #doctest: +ELLIPSIS
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
            featuretypes.update(i[2])
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

        for feature, featuretypes in self.feature_classes.items():
            # get rid of "unannotated"
            ft = featuretypes.difference(['.'])
            key = frozenset(ft)
            self.class_features[key].append(feature)
            self.class_counts[key] += 1

        # convert defaultdicts to regular dicts
        self.class_features = dict(self.class_features)
        self.class_counts = dict(self.class_counts)
        self.feature_classes = dict(self.feature_classes)
