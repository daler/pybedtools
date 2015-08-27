_bedtools_path = ""
_tabix_path = ""
_bgzip_path = ""
_R_path = ""

tempfile_prefix = 'pybedtools.'
tempfile_suffix = '.tmp'


# Checking for BEDTools will happen when creating the first BedTool; other
# checks happen at first use (BAM object creation; tabix-ing a BedTool)
_bedtools_installed = False
_tabix_installed = False
_bgzip_installed = False
_R_installed = False
_v_2_15_plus = False

KEEP_TEMPFILES = False
_DEBUG = True

# Check calls against these names to only allow calls to known BEDTools
# programs (basic security)
#
_prog_names = {

    # Genome arithmetic
    'intersectBed': 'intersect',
    'windowBed': 'window',
    'closestBed': 'closest',
    'coverageBed': 'coverage',
    'mapBed': 'map',
    'genomeCoverageBed': 'genomecov',
    'mergeBed': 'merge',
    'clusterBed': 'cluster',
    'complementBed': 'complement',
    'subtractBed': 'subtract',
    'slopBed': 'slop',
    'flankBed': 'flank',
    'sortBed': 'sort',
    'randomBed': 'random',
    'shuffleBed': 'shuffle',
    'annotateBed': 'annotate',

    # multi-way
    'multiIntersectBed': 'multiinter',
    'unionBedGraphs': 'unionbedg',

    # PE
    'pairToBed': 'pairtobed',
    'pairToPair': 'pairtopair',

    # format conversion
    'bamToBed': 'bamtobed',
    'bedToBam': 'bedtobam',
    'bedpeToBam': 'bedpetobam',
    'bed12ToBed6': 'bed12tobed6',
    'bamToFastq': 'bamtofastq',

    # fasta
    'fastaFromBed': 'getfasta',
    'maskFastaFromBed': 'maskfasta',
    'nucBed': 'nuc',

    # bam-centric
    'multiBamCov': 'multicov',
    'tagBam': 'tag',

    # stats
    'jaccard': 'jaccard',
    'reldist': 'reldist',

    # misc
    'getOverlap': 'overlap',
    'bedToIgv': 'igv',
    'linksBed': 'links',
    'windowMaker': 'makewindows',
    'groupBy': 'groupby',
    'expandCols': 'expand',
    'sample': 'sample',
    'fisher': 'fisher',
}

_old_names = list(_prog_names.keys())
_new_names = list(_prog_names.values())

_column_names = {
    'bed': ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart',
            'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'],
    'gff': ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'],
}
