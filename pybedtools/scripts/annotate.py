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
from pybedtools.cbedtools import parse_attributes, create_interval_from_list
import collections

def get_gff_name(field):
    attrs = parse_attributes(field)
    for key in ("ID", "gene_name", "transcript_id", "gene_id", "Parent"):
        if key in attrs: return attrs[key]

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
    a, full_b = BedTool(aname), BedTool(bname)

    def exonify(f):
        f[4] = "exon"
        return f

    bintrons = full_b.introns()
    bexons = full_b.bed6().each(exonify)

    fb = open(BedTool._tmp(), "w")
    for bi in bintrons:
        fb.write(str(bi) + "\n")
    for be in bexons:
        fb.write(str(be) + "\n")
    fb.close()

    b = BedTool(fb.name).sort()
    b.saveas('b.bed')
    afields = a.field_count()
    c = a.closest(b, d=True, t="all")
    get_name = gen_get_name(b, afields)

    dbed = open(BedTool._tmp(), "w")

    # keep the name, distance and feature type.
    seen_by_line = collections.defaultdict(list)
    c.saveas('c.bed')
    assert len(c) > 1
    for feat in c:
        fields = feat.fields
        key = "\t".join(fields[:afields])
        distance = int(fields[-1])
        seen_by_line[key].append({'distance': distance, 'name': get_name(feat),
                                  'type': fields[afields + 4]})


    for key, dist_info in seen_by_line.iteritems():
        if len(dist_info) > 1:
            assert all(d['distance'] == 0 for d in dist_info), (dist_info)

        # if the distance is zero. figure it it's intron or exon.
        dist = dist_info[0]['distance']
        if all(d['distance'] == 0 for d in dist_info):
            dist = ";".join(sorted(set(d['type'] for d in dist_info)))

        names = ",".join(sorted(set(d['name'] for d in dist_info)))
        new_line = "\t".join([key, names, str(dist)])
        dbed.write(new_line + "\n")
    dbed.close()
    d = BedTool(dbed.name)
    d.saveas('y.bed')
    assert len(d) == len(a)
    return d

def add_xstream(a, b, dist, updown, report_distance=False):
    # run a window up or downstream.
    dir = dict(up="l", down="r")[updown]
    kwargs = {'sw':True, dir: dist}

    # have to set the other to 0
    if "l" in kwargs: kwargs["r"] = 0
    else: kwargs["l"] = 0

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
        if key in seen: continue
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
