"""
annotate a file with the neearest features in another.

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
from pybedtools.cbedtools import Attributes
import collections

# PYTHONPATH=$PYTHONPATH:. python scripts/annotate.py -a data/new.regions.bed
# -b data/Homo_sapiens.hg18.gtf --upstream 5000

# $ pybedtools annotate -a regions.bed -b knownGene.bed --upstream 10000
#                  --downstream 5000 --report-distance
# a bed: regions.bed and another:
# annotation.bed, it would add 4 columns to regions.bed:
# nearest-feature, nearest-distance, upstream-features, downstream-features
# where the up/downstream features are determined by a distance
# parameter, e.g. like --upstream 10000 --downstream 5000


def get_gff_name(field):
    attrs = Attributes(field)
    for key in ("ID", "gene_name", "transcript_id", "gene_id", "Parent"):
        if key in attrs:
            return attrs[key]


def gen_get_name(b, afields):
    btype = b.file_type
    if btype == "bed":
        get_name = lambda fields: fields[afields + 3]
    elif btype == "gff":
        def get_name(fields):
            return get_gff_name(fields[afields + 7])
    else:
        raise Exception("not implemented")
    return get_name


def add_closest(aname, bname):
    a, b = BedTool(aname), BedTool(bname)

    afields = a.field_count()
    c = a.closest(b, d=True)
    get_name = gen_get_name(b, afields)

    dbed = open(BedTool._tmp(), "w")
    # keep the name and distance
    seen_by_line = collections.defaultdict(list)
    for feat in c:
        key = "\t".join(feat[:afields])
        seen_by_line[key].append([feat[-1], get_name(feat)])

    for key, dist_names in seen_by_line.items():
        if len(dist_names) > 0:
            assert len(set([d[0] for d in dist_names])) == 1
        names = ",".join(sorted(set(d[1] for d in dist_names)))
        new_line = "\t".join([key] + [names] + [dist_names[0][0]])
        dbed.write(new_line + "\n")
    dbed.close()
    d = BedTool(dbed.name)
    assert len(d) == len(a)
    return d


def add_xstream(a, b, dist, updown, report_distance=False):
    # run a window up or downstream.
    dir = dict(up="l", down="r")[updown]
    kwargs = {'sw': True, dir: dist}

    # have to set the other to 0
    if "l" in kwargs:
        kwargs["r"] = 0
    else:
        kwargs["l"] = 0

    c = a.window(b, **kwargs)
    afields = a.field_count()

    get_name = gen_get_name(b, afields)

    seen = collections.defaultdict(set)
    # condense to unique names.
    for row in c:
        key = "\t".join(row[:afields])
        seen[key].update([get_name(row)])

    d = open(BedTool._tmp(), "w")
    for row in seen:
        d.write(row + "\t" + ",".join(sorted(seen[row])) + "\n")

    # write the entries that did not appear in the window'ed Bed
    for row in a:
        key = "\t".join(row[:afields])
        if key in seen:
            continue
        d.write(str(row) + "\t.\n")

    d.close()
    dbed = BedTool(d.name)
    assert len(dbed) == len(a)
    return dbed


def main():
    """
    annotate a file with the neearest features in another.
    """
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
    b = BedTool(args.b)
    # TODO: support --report-distance for up/downstream.
    if args.upstream:
        c = add_xstream(c, b, args.upstream, "up", args.report_distance)
    if args.downstream:
        c = add_xstream(c, b, args.downstream, "down",
                        args.report_distance)

    for row in c.sort():
        print(row)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS).failed == 0:
        main()
