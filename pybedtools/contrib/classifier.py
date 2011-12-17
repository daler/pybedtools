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

        d = defaultdict(set)

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

            key = '\t'.join(feature[:bed_fields])
            d[key].update([featuretype])

        return d

    def class_counts(self, d):
        """
        `d` is the dictionary returned by classify.  Returns a dictionary of
        featuretype classes as keys and the number of BED features in the class
        as values.

        A featuretype class can be one featuretype or a combination (e.g.,
        "exon and intron")
        """
        count_d = defaultdict(int)
        for feature, featuretypes in d.items():
            # get rid of "unannotated"
            ft = featuretypes.difference(['.'])
            label = ','.join(sorted(list(ft)))
            count_d[label] += 1

        return count_d
