"""
Chromsize dictionaries, as downloaded from UCSC.  Need more? Use::

    pybedtools.get_chromsizes_from_ucsc('assemblyname')

"""
# Figure out which version of OrderedDict we want....
import sys
if (sys.version_info[0] == 2) and (sys.version_info[1] < 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict


dm3 = OrderedDict((
    ('chr2L', (0, 23011544)),
    ('chr2R', (0, 21146708)),
    ('chr3L', (0, 24543557)),
    ('chr3R', (0, 27905053)),
    ('chr4', (0, 1351857)),
    ('chrX', (0, 22422827)),
    ('chr2LHet', (0, 368872)),
    ('chr2RHet', (0, 3288761)),
    ('chr3LHet', (0, 2555491)),
    ('chr3RHet', (0, 2517507)),
    ('chrM', (0, 19517)),
    ('chrU', (0, 10049037)),
    ('chrUextra', (0, 29004656)),
    ('chrXHet', (0, 204112)),
    ('chrYHet', (0, 347038)),
))

# No chrUextra or chrM
dm3.default = OrderedDict()
for chrom, size in dm3.items():
    if chrom in ['chrUextra', 'chrM']:
        continue
    dm3.default[chrom] = size

# No chrU*, chr*Het, or chrM
dm3.euchromatic = OrderedDict()
for chrom, size in dm3.default.items():
    if 'chrU' in chrom:
        continue
    if 'Het' in chrom:
        continue
    dm3.euchromatic[chrom] = size


mm9 = OrderedDict((
    ('chr1', (0, 197195432)),
    ('chr2', (0, 181748087)),
    ('chr3', (0, 159599783)),
    ('chr4', (0, 155630120)),
    ('chr5', (0, 152537259)),
    ('chr6', (0, 149517037)),
    ('chr7', (0, 152524553)),
    ('chr8', (0, 131738871)),
    ('chr9', (0, 124076172)),
    ('chr10', (0, 129993255)),
    ('chr11', (0, 121843856)),
    ('chr12', (0, 121257530)),
    ('chr13', (0, 120284312)),
    ('chr14', (0, 125194864)),
    ('chr15', (0, 103494974)),
    ('chr16', (0, 98319150)),
    ('chr17', (0, 95272651)),
    ('chr18', (0, 90772031)),
    ('chr19', (0, 61342430)),
    ('chrX', (0, 166650296)),
    ('chrY', (0, 15902555)),
    ('chrM', (0, 16299)),
    ('chr13_random', (0, 400311)),
    ('chr16_random', (0, 3994)),
    ('chr17_random', (0, 628739)),
    ('chr1_random', (0, 1231697)),
    ('chr3_random', (0, 41899)),
    ('chr4_random', (0, 160594)),
    ('chr5_random', (0, 357350)),
    ('chr7_random', (0, 362490)),
    ('chr8_random', (0, 849593)),
    ('chr9_random', (0, 449403)),
    ('chrUn_random', (0, 5900358)),
    ('chrX_random', (0, 1785075)),
    ('chrY_random', (0, 58682461)),
))

mm9.default = OrderedDict()
for chrom, size in mm9.items():
    if '_random' in chrom:
        continue
    if chrom == 'chrM':
        continue
    mm9.default[chrom] = size


hg18 = OrderedDict((
    ('chr1', (0, 247249719)),
    ('chr2', (0, 242951149)),
    ('chr3', (0, 199501827)),
    ('chr4', (0, 191273063)),
    ('chr5', (0, 180857866)),
    ('chr6', (0, 170899992)),
    ('chr7', (0, 158821424)),
    ('chr8', (0, 146274826)),
    ('chr9', (0, 140273252)),
    ('chr10', (0, 135374737)),
    ('chr11', (0, 134452384)),
    ('chr12', (0, 132349534)),
    ('chr13', (0, 114142980)),
    ('chr14', (0, 106368585)),
    ('chr15', (0, 100338915)),
    ('chr16', (0, 88827254)),
    ('chr17', (0, 78774742)),
    ('chr18', (0, 76117153)),
    ('chr19', (0, 63811651)),
    ('chr20', (0, 62435964)),
    ('chr21', (0, 46944323)),
    ('chr22', (0, 49691432)),
    ('chrX', (0, 154913754)),
    ('chrY', (0, 57772954)),
    ('chrM', (0, 16571)),
    ('chr10_random', (0, 113275)),
    ('chr11_random', (0, 215294)),
    ('chr13_random', (0, 186858)),
    ('chr15_random', (0, 784346)),
    ('chr16_random', (0, 105485)),
    ('chr17_random', (0, 2617613)),
    ('chr18_random', (0, 4262)),
    ('chr19_random', (0, 301858)),
    ('chr1_random', (0, 1663265)),
    ('chr21_random', (0, 1679693)),
    ('chr22_h2_hap1', (0, 63661)),
    ('chr22_random', (0, 257318)),
    ('chr2_random', (0, 185571)),
    ('chr3_random', (0, 749256)),
    ('chr4_random', (0, 842648)),
    ('chr5_h2_hap1', (0, 1794870)),
    ('chr5_random', (0, 143687)),
    ('chr6_cox_hap1', (0, 4731698)),
    ('chr6_qbl_hap2', (0, 4565931)),
    ('chr6_random', (0, 1875562)),
    ('chr7_random', (0, 549659)),
    ('chr8_random', (0, 943810)),
    ('chr9_random', (0, 1146434)),
    ('chrX_random', (0, 1719168)),
))

hg18.default = OrderedDict()
for chrom, size in hg18.items():
    if '_' in chrom:
        continue
    if chrom == 'chrM':
        continue
    hg18.default[chrom] = size


hg19 = OrderedDict((
    ('chr1', (0, 249250621)),
    ('chr2', (0, 243199373)),
    ('chr3', (0, 198022430)),
    ('chr4', (0, 191154276)),
    ('chr5', (0, 180915260)),
    ('chr6', (0, 171115067)),
    ('chr7', (0, 159138663)),
    ('chr8', (0, 146364022)),
    ('chr9', (0, 141213431)),
    ('chr10', (0, 135534747)),
    ('chr11', (0, 135006516)),
    ('chr12', (0, 133851895)),
    ('chr13', (0, 115169878)),
    ('chr14', (0, 107349540)),
    ('chr15', (0, 102531392)),
    ('chr16', (0, 90354753)),
    ('chr17', (0, 81195210)),
    ('chr18', (0, 78077248)),
    ('chr19', (0, 59128983)),
    ('chr20', (0, 63025520)),
    ('chr21', (0, 48129895)),
    ('chr22', (0, 51304566)),
    ('chrX', (0, 155270560)),
    ('chrY', (0, 59373566)),
    ('chrM', (0, 16571)),
    ('chr6_ssto_hap7', (0, 4928567)),
    ('chr6_mcf_hap5', (0, 4833398)),
    ('chr6_cox_hap2', (0, 4795371)),
    ('chr6_mann_hap4', (0, 4683263)),
    ('chr6_apd_hap1', (0, 4622290)),
    ('chr6_qbl_hap6', (0, 4611984)),
    ('chr6_dbb_hap3', (0, 4610396)),
    ('chr17_ctg5_hap1', (0, 1680828)),
    ('chr4_ctg9_hap1', (0, 590426)),
    ('chr1_gl000192_random', (0, 547496)),
    ('chrUn_gl000225', (0, 211173)),
    ('chr4_gl000194_random', (0, 191469)),
    ('chr4_gl000193_random', (0, 189789)),
    ('chr9_gl000200_random', (0, 187035)),
    ('chrUn_gl000222', (0, 186861)),
    ('chrUn_gl000212', (0, 186858)),
    ('chr7_gl000195_random', (0, 182896)),
    ('chrUn_gl000223', (0, 180455)),
    ('chrUn_gl000224', (0, 179693)),
    ('chrUn_gl000219', (0, 179198)),
    ('chr17_gl000205_random', (0, 174588)),
    ('chrUn_gl000215', (0, 172545)),
    ('chrUn_gl000216', (0, 172294)),
    ('chrUn_gl000217', (0, 172149)),
    ('chr9_gl000199_random', (0, 169874)),
    ('chrUn_gl000211', (0, 166566)),
    ('chrUn_gl000213', (0, 164239)),
    ('chrUn_gl000220', (0, 161802)),
    ('chrUn_gl000218', (0, 161147)),
    ('chr19_gl000209_random', (0, 159169)),
    ('chrUn_gl000221', (0, 155397)),
    ('chrUn_gl000214', (0, 137718)),
    ('chrUn_gl000228', (0, 129120)),
    ('chrUn_gl000227', (0, 128374)),
    ('chr1_gl000191_random', (0, 106433)),
    ('chr19_gl000208_random', (0, 92689)),
    ('chr9_gl000198_random', (0, 90085)),
    ('chr17_gl000204_random', (0, 81310)),
    ('chrUn_gl000233', (0, 45941)),
    ('chrUn_gl000237', (0, 45867)),
    ('chrUn_gl000230', (0, 43691)),
    ('chrUn_gl000242', (0, 43523)),
    ('chrUn_gl000243', (0, 43341)),
    ('chrUn_gl000241', (0, 42152)),
    ('chrUn_gl000236', (0, 41934)),
    ('chrUn_gl000240', (0, 41933)),
    ('chr17_gl000206_random', (0, 41001)),
    ('chrUn_gl000232', (0, 40652)),
    ('chrUn_gl000234', (0, 40531)),
    ('chr11_gl000202_random', (0, 40103)),
    ('chrUn_gl000238', (0, 39939)),
    ('chrUn_gl000244', (0, 39929)),
    ('chrUn_gl000248', (0, 39786)),
    ('chr8_gl000196_random', (0, 38914)),
    ('chrUn_gl000249', (0, 38502)),
    ('chrUn_gl000246', (0, 38154)),
    ('chr17_gl000203_random', (0, 37498)),
    ('chr8_gl000197_random', (0, 37175)),
    ('chrUn_gl000245', (0, 36651)),
    ('chrUn_gl000247', (0, 36422)),
    ('chr9_gl000201_random', (0, 36148)),
    ('chrUn_gl000235', (0, 34474)),
    ('chrUn_gl000239', (0, 33824)),
    ('chr21_gl000210_random', (0, 27682)),
    ('chrUn_gl000231', (0, 27386)),
    ('chrUn_gl000229', (0, 19913)),
    ('chrUn_gl000226', (0, 15008)),
    ('chr18_gl000207_random', (0, 4262)),
))

hg19.default = OrderedDict()
for chrom, size in hg19.items():
    if '_' in chrom:
        continue
    if chrom == 'chrM':
        continue
    hg19.default[chrom] = size
