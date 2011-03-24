
class Feature(object):
    line_sep = "\t"
    position_locs = (0, 1, 2)
    __slots__ = ('chr', 'start', 'stop', '_line_arr')
    def __init__(self, line):

        line_arr = line.rstrip("\r\n").split(self.line_sep)
        locs = self.__class__.position_locs
        self.chr = line_arr[locs[0]]
        self.start = int(line_arr[locs[1]])
        self.stop = int(line_arr[locs[2]])
        self._line_arr = line_arr

        self._post_initialize()

    def __getattr__(self, field):
        try:
            idx = self.__slots__.index(field)
        except ValueError:
            raise AttributeError("%s does not exist" % field)
        return self._line_arr[idx]

    def _post_initialize(self):
        pass

    def __repr__(self):
        return "%s(%s:%i-%i)" % (self.__class__.__name__, self.chr, self.start, self.stop)

class BedFeature(Feature):
    r"""
    >>> b = BedFeature("chr1\t23\t45\thello")
    >>> b.name
    'hello'
    >>> b = BedFeature("chr10\t99\t122\tgoodbye")
    >>> b
    BedFeature(chr10:99-122)
    """
    __slots__ = ('chr', 'start', 'stop', 'name', 'value', 'strand',
                 'thickStart', 'thickStop', 'itemRGB',
                 'blockCount', 'blockSizes', 'blockStarts')

class GFFFeature(Feature):
    position_locs = (0, 3, 4)
    __slots__ = ('chr', 'method', 'featuretype', 'start', 'stop', 'score',
                 'strand', 'phase', '_sattributes', '_attributes')
    _sep = ";"
    _field_sep = "="

    @property
    def attributes(self):
        return self._attributes

    def _post_initialize(self):
        # copied from genomicfeatures
        # this could also be done lazily the first time attributes() is called.
        # i think that's a better option in the spirit of keeping it minimal.
        _attributes = {}
        kvs = map(str.strip, self._sattributes.strip().split(self._sep))
        for field, value in [kv.split(self._field_sep) for kv in kvs if kv]:
            _attributes[field] = value.replace('"', '')
        self._attributes = _attributes

    @property
    def name(self):
        return self.attributes.get("ID", self.attributes.get("gene_name"))


class GTFFeature(GFFFeature):
    r"""
    >>> g = GTFFeature('11\tpseudogene\texon\t75780\t76143\t.\t+\t.\tgene_id "ENSG00000253826"; transcript_id "ENST00000519787"; exon_number "1"; gene_name "RP11-304M2.1"; transcript_name "RP11-304M2.1-201";')
    >>> g
    GTFFeature(11:75780-76143)
    >>> g.name
    'RP11-304M2.1'
    >>> g.strand
    '+'

    >>> g = GTFFeature('11\tpseudogene\texon\t75780\t76143\t.\t-\t.\t')
    >>> g.strand, g.score, g.phase
    ('-', '.', '.')
    """
    _sep = ';'
    _field_sep = " "
    __slots__ = GFFFeature.__slots__


if __name__ == "__main__":
    import doctest
    doctest.testmod()
