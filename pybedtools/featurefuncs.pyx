# distutils: language = c++
from cbedtools cimport Interval
from cbedtools import create_interval_from_list


cpdef extend_fields(Interval feature, int n):
    """
    Pads the fields of the feature with "." to a total length of `n` fields,
    """
    fields = feature.fields[:]
    while len(fields) < n:
        fields.append('.')
    i = create_interval_from_list(fields)

    if n > 4 and (i[4] == '.'):
        i[4] = '0'
    if n > 6 and (i[6] == '.'):
        i[6] = str(i.start)
    if n > 7 and (i[7] == '.'):
        i[7] = str(i.stop)
    if n > 8 and (i[8] == '.'):
        i[8] = '0,0,0'
    return i



cpdef center(Interval feature, int width=100):
    """
    Return the *width* bp from the center of a feature.  If a feature is
    smaller than *width*, then return the entire feature.
    """
    if len(feature) < width:
        return feature
    cdef int start = feature.start
    cdef int stop = feature.stop
    cdef int center = start + (stop - start) / 2
    halfwidth = width / 2
    feature.start = center - halfwidth
    if feature.start < 1:
        feature.start = 1
    if halfwidth == 0:
        halfwidth = 1
    feature.stop = center + halfwidth
    return feature


cpdef midpoint(Interval feature):
    """
    Specialized version of `center()` that just returns the single-bp midpoint
    """
    start = feature.start + (feature.stop - feature.start) / 2
    stop = start + 1
    feature.start = start
    feature.stop = stop
    return feature


cpdef greater_than(Interval feature, int size=100):
    """
    Return True if feature length > *size*
    """
    return len(feature) > size


cpdef less_than(Interval feature, int size=100):
    """
    Return True if feature length < *size*
    """
    return len(feature) < size


cpdef normalized_to_length(Interval feature, int idx=4, float scalar=0.001):
    """
    Normalizes the value at feature[idx] to the feature's length, in kb.

    *idx*, by default, is the score field for a BED file, but specify any
    integer.

    The value at *idx* will be replaced with its scaled value.

    *scalar* will be multiplied by the value at *idx*, by default this is
    0.001, or per kb.

    Useful for calculating RPKM after running intersect with counts
    """
    feature[idx] = str(float(feature[idx]) * scalar / len(feature))
    return feature


cpdef rename(Interval feature, str name):
    """
    Forces a rename of all features, e.g., for renaming everything in a file
    'exon'
    """
    feature.name = name
    return feature


cpdef bedgraph_scale(Interval feature, float scalar):
    feature[3] = str(float(feature[3]) * scalar)
    return feature


cpdef TSS(Interval feature, int upstream=500, int downstream=500, add_to_name=None, genome=None):
    """
    Alias for five_prime.
    """
    return star_prime(feature, upstream, downstream, prime=5,
                      add_to_name=add_to_name, genome=genome)


cdef star_prime(Interval feature, int upstream=500, int downstream=500, int prime=5,
                 add_to_name=None, genome=None):

    if prime == 5:
        if feature.strand == '-':
            start = feature.stop - downstream
            stop = feature.stop + upstream
        else:
            start = feature.start - upstream
            stop = feature.start + downstream
    elif prime == 3:
        if feature.strand == '-':
            start = feature.start - downstream
            stop = feature.start + upstream
        else:
            start = feature.stop - upstream
            stop = feature.stop + downstream
    if add_to_name:
        try:
            feature.name += add_to_name
        except AttributeError:
            pass
    if genome is not None:
        gstart, gstop = genome[feature.chrom]
        stop = min(stop, gstop)
        start = max(start, gstart)
    if start < 0:
        start = 0
    if start > stop:
        start = stop
    feature.start = start
    feature.stop = stop
    return feature

cpdef five_prime(Interval feature, int upstream=500, int downstream=500,
                 add_to_name=None, genome=None):
    """
    Returns the 5'-most coordinate, plus `upstream` and `downstream` bp; adds
    the string `add_to_name` to the feature's name if provided (e.g., "_TSS")

    Parameters
    ----------
    feature : pybedtools.Interval instance

    upstream, downstream : int
        Number of bp upstream or downstream of the strand-specific start
        position of the feature to include. Default is 500 for both upstream
        and downstream so that the returned feature is 1kb centered on the 5'
        end of the feature. Unstranded features (where strand=".") are treated
        as plus-strand features.

    add_to_name : str or None
        If not None, append the string suffix to the name field of the feature (for
        example "_TSS").

    genome : dict or None
        If not None, then ensure that the start/stop positions are within the
        boundaries of the chromosome.
    """
    return star_prime(feature, upstream, downstream, prime=5,
                      add_to_name=add_to_name, genome=genome)


cpdef three_prime(Interval feature, int upstream=500, int downstream=500,
                  add_to_name=None, genome=None):
    """
    Returns the 3'-most coordinate, plus `upstream` and `downstream` bp; adds
    the string `add_to_name` to the feature's name if provided (e.g.,
    "_polyA_site")

    Parameters
    ----------
    feature : pybedtools.Interval instance

    upstream, downstrea : int
        Number of bp upstream or downstream of the strand-specific stop
        position of the feature to include. Default is 500 for both upstream
        and downstream so that the returned feature is 1kb centered on the 5'
        end of the feature. Unstranded features (where strand=".") are treated
        as plus-strand features.

    add_to_name : str or None
        If not None, append the string suffix to the name field of the feature (for
        example "_TSS").

    genome : dict or None
        If not None, then ensure that the start/stop positions are within the
        boundaries of the chromosome.


    """
    return star_prime(feature, upstream, downstream, prime=3,
                      add_to_name=add_to_name, genome=genome)

cpdef add_color(Interval feature, cmap, norm):
    """
    Signature:

        add_color(feature, cmap, norm)

    Given the matplotlib colormap `cmap` and the matplotlib Normalize instance
    `norm`, return a new 9-field feature (extended out if needed) with the RGB
    tuple set according to the score.
    """
    if len(feature.fields) < 9:
        feature = extend_fields(feature, 9)
        feature[6] = str(feature.start)
        feature[7] = str(feature.stop)

    rgb_float = cmap(norm(float(feature.score)))
    feature[8] = ','.join([str(int(i * 255)) for i in rgb_float[:3]])
    return feature


cpdef gff2bed(Interval feature, name_field=None):
    """
    Signature:

        gff2bed(feature, name_field=None)

    Converts a GFF feature into a BED6 feature.  By default, the name of the
    new BED will be feature.name, but if `name_field` is provided then the name
    of the new BED will be feature.attrs[name_field].

    `name_field` can also be an integer to index into the fields of the object,
    so if you want the BED name to be the GFF featuretype, then you can use
    `name_field=2`.

    If the specified field does not exist, then "." will be used for the name.
    """
    if name_field is None:
        name = feature.name
    else:
        try:
            if isinstance(name_field, basestring):
                name = feature.attrs[name_field]
            if isinstance(name_field, int):
                name = feature[name_field]
        except (NameError, KeyError):
            name = "."
    return create_interval_from_list([
        str(feature.chrom),
        str(feature.start),
        str(feature.stop),
        name,
        feature.score,
        feature.strand])


cpdef bed2gff(Interval feature):
    """
    Signature:

        bed2gff(feature)

    Converts a BED feature (BED3 through BED12) into a GFF format.

    Chrom, start, stop, score, and strand are put directly into the
    corresponding GFF fields.  Other BED fields are put into the GFF attributes
    field, named according to the UCSC BED format definition.

    If there are more than 12 BED fields, the additional fields will be added
    to the GFF attributes using the 0-based index (so starting at "12") as the
    key.

    GFF fields that do not have a direct mapping to BED format (feature type,
    source, phase) are set to ".".

    1 bp is added to the start position to finish the conversion to GFF.
    """

    # Note that Interval.score, .strand, and .name have a default of ".", so no
    # need to do the extra try/except IndexError for those fields.
    mapping = (
        (6, "thickStart"),
        (7, "thickEnd"),
        (8, "itemRgb"),
        (9, "blockCount"),
        (10, "blockSizes"),
        (11, "blockStarts")
    )

    # Add any standard BED fields we might have
    attributes = ['Name="%s"' % feature.name]
    for k, v in mapping:
        try:
            attributes.append('%s="%s"' % (v, feature.fields[k]))
        except IndexError:
            break

    # Add any additional fields, keyed by their index
    if len(feature.fields) > 12:
        for i in range(12, len(feature.fields)):
            attributes.append('%s="%s"' % (i, feature.fields[i]))

    attributes = '; '.join(attributes) + ';'

    return create_interval_from_list([
        str(feature.chrom),
        '.',
        '.',
        str(feature.start + 1),
        str(feature.stop),
        feature.score,
        feature.strand,
        '.',
        attributes])


class UniqueID(object):
    def __init__(self, pattern="%d", first=0):
        """
        Class to help create uniquely-named features.

        Example usage:

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> uid = UniqueID("f_%d")
        >>> print(a.each(uid))  # doctest: +NORMALIZE_WHITESPACE
        chr1    1    100    f_0    0    +
        chr1  100    200    f_1    0    +
        chr1  150    500    f_2    0    -
        chr1  900    950    f_3    0    +

        Parameters
        ----------
        pattern : str

            Pattern will be filled in using `% self.count`

        first : int
            `self.count` will be initialzed to this value.

        """
        self.pattern = pattern
        self.first = first
        self.count = first

    def __call__(self, feature):
        feature.name = self.pattern % self.count
        self.count += 1
        return feature

