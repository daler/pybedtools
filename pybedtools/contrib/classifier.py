import pybedtools
import itertools
from collections import defaultdict


class Classifier(object):
    def __init__(self, bed, annotations):
        """
        Class for working with peaks in overlapping annotated regions --
        introns, exons, etc
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
        Intersect the BED file with the annotations file and return
        a dictionary where keys are BED features and values are the set of
        featuretypes that BED feature was found in.


        * `include` is an optional list of featuretypes to restrict the
          classification to
        * `exclude` is an optional list of featuretypes to
          exclude from classification (all other featuretypes will be used).

        To see what's available, use available_featuretypes().

        When run, this method creates the following dictionaries as attributes
        of this object:

         feature_classes: keys are Intervals from `bed`;
                          values are sets of featuretypes from `annotations`

         class_features : keys are frozensets of featuretypes from
                          `annotations`; values are lists of Intervals from
                          `bed`;

         class_counts   : keys are frozensets of featuretypes from
                          annotations`; values are number of features -- so the
                          length of values in the class_features dictionary.

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
