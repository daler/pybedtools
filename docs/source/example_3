import sys
import multiprocessing
import pybedtools

# get example GFF and BAM filenames
gff = pybedtools.example_filename('gdc.gff')
bam = pybedtools.example_filename('gdc.bam')

# Some GFF files have invalid entries -- like chromosomes with negative coords
# or features of length = 0.  This line removes them and saves the result in a
# tempfile
g = pybedtools.BedTool(gff).remove_invalid().saveas()


# Next, we create a function to pass only features for a particular
# featuretype.  This is similar to a "grep" operation when applied to every
# feature in a BedTool
def featuretype_filter(feature, featuretype):
    if feature[2] == featuretype:
        return True
    return False


# This function will eventually be run in parallel, applying the filter above
# to several different BedTools simultaneously
def subset_featuretypes(featuretype):
    result = g.filter(featuretype_filter, featuretype).saveas()
    return pybedtools.BedTool(result.fn)


# This function performs the intersection of a BAM file with a GFF file and
# returns the total number of hits.  It will eventually be run in parallel.
def count_reads_in_features(features_fn):
    """
    Callback function to count reads in features
    """
    # BAM files are auto-detected; no need for an `abam` argument.  Here we
    # construct a new BedTool out of the BAM file and intersect it with the
    # features filename.

    # We use stream=True so that no intermediate tempfile is
    # created, and bed=True so that the .count() method can iterate through the
    # resulting streamed BedTool.
    return pybedtools.BedTool(bam).intersect(
                             b=features_fn,
                             stream=True).count()


# Set up a pool of workers for parallel processing
pool = multiprocessing.Pool()

# Create separate files for introns and exons, using the function we defined
# above
featuretypes = ('intron', 'exon')
introns, exons = pool.map(subset_featuretypes, featuretypes)

# Perform some genome algebra to get unique and shared intron/exon regions.
# Here we keep only the filename of the results, which is safer than an entire
# BedTool for passing around in parallel computations.
exon_only = exons.subtract(introns).merge().remove_invalid().saveas().fn
intron_only = introns.subtract(exons).merge().remove_invalid().saveas().fn
intron_and_exon = exons.intersect(introns).merge().remove_invalid().saveas().fn

# Do intersections with BAM file in parallel, using the other function we
# defined above
features = (exon_only, intron_only, intron_and_exon)
results = pool.map(count_reads_in_features, features)

# Print the results
labels = ('      exon only:',
          '    intron only:',
          'intron and exon:')

for label, reads in zip(labels, results):
    sys.stdout.write('%s %s\n' % (label, reads))
