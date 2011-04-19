"""
    %prog [options]

if --upstream and --downstream are not specified only 2 colummns are
added to the [a] file: nearest_name, nearest_distance.
if they are specified upstream_names, downstream_names are specified.
if --report-distance is included, the format will be:

   gene:dist,gene2:dist2

for each of the up and downstream columns.
"""

import argparse
import sys
from pybedtools import BedTool
from pybedtools.helpers import parse_attributes

# $ pybedtools annotate -a regions.bed -b knownGene.bed --upstream 10000
#                  --downstream 5000 --report-distance
# a bed: regions.bed and another:
# annotation.bed, it would add 4 columns to regions.bed:
# nearest-feature, nearest-distance, upstream-features, downstream-features
# where the up/downstream features are determined by a distance
# parameter, e.g. like --upstream 10000 --downstream 5000

def add_closest(aname, bname):
    a, b = BedTool(aname), BedTool(bname)

    afields = a.field_count()
    c = a.closest(b, d=True)
    btype = b.file_type
    if btype == "bed":
        get_name = lambda fields: fields[afields + 4]
    elif btype == "gff":
        def get_name(fields):
            attrs = parse_attributes(fields[afields + 8])
            for key in ("ID", "gene_name", "transcript_id", "gene_id",
                                                            "Parent"):
                if key in attrs: return attrs[key]
    else:
        raise Exception("not implemented")

    d = open(c._tmp(), "w")
    # keep the name and distance
    seen = {}
    for row in c:
        new_line = "\t".join(row[:afields] + [get_name(row), row[-1]])
        if new_line in seen: continue
        print >>d, new_line
        seen[new_line] = None
    d.close()
    return BedTool(d.name)


def main():
    p = argparse.ArgumentParser(description=__doc__, prog=sys.argv[0])
    p.add_argument("-a", dest="a", help="file to annotate")
    p.add_argument("-b", dest="b", help="file with annotations")
    p.add_argument("--upstream", dest="upstream", type=int, default=None,
                   help="distance upstream of [a] to look for [b]")
    p.add_argument("--downstream", dest="downstream", type=int, default=None,
                   help="distance downstream of [a] to look for [b]")
    p.add_argument("--report-distance", dest="report_distance", default=False,
                   help="report the distance, not just the genes",
                   action="store_true")
    args = p.parse_args()
    if (args.a is None or args.b is None):
        sys.exit(not p.print_help())

    c = add_closest(args.a, args.b)
    c.saveas('g.bed')

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS).failed == 0:
        main()
