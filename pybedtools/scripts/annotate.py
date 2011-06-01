#!/usr/bin/env python
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
    """
    Returns the name, parsed from the GFF attributes field, *field*
    """
    attrs = parse_attributes(field)
    for key in ("ID", "gene_name", "transcript_id", "gene_id", "Parent"):
        if key in attrs:
            return attrs[key]


def gen_get_name(b, afields):
    """
    Factory for making name-getter functions in a format-specific way
    """
    btype = b.file_type
    if btype == "bed":
        get_name = lambda fields: fields[afields + 3]
    elif btype == "gff":
        def get_name(fields):
            return get_gff_name(fields[afields + 7])
    else:
        raise Exception("not implemented")
    return get_name


def add_closest(aname, bname, save_intermediates=False):
    """
    Tags each feature in *aname* with the closest feature, or if the closest
    feature had a distance of 0, figure out if it was in an intron or exon.
    """
    a, full_b = BedTool(aname), BedTool(bname)

    def exonify(f):
        f[4] = "exon"
        return f

    # Construct introns from a BED12
    bintrons = full_b.introns()

    # From that same BED12, split out the exons and put 'exon' in the score
    # field
    bexons = full_b.bed6().each(exonify)

    fb = open(BedTool._tmp(), "w")
    for bi in bintrons:
        fb.write(str(bi) + "\n")
    for be in bexons:
        fb.write(str(be) + "\n")
    fb.close()

    if save_intermediates:
        bfn = 'b.bed'
    else:
        bfn = None
    b = BedTool(fb.name).sort().saveas(bfn)

    # Alternatively:
    # b = bintrons.cat(bexons).sort().saveas()

    afields = a.field_count()

    # d=True will report distance to closest feature
    # t='all' will report all ties instead of choosing just one to report.
    if save_intermediates:
        cfn = 'c.bed'
    else:
        cfn = None
    c = a.closest(b, d=True, t="all").saveas(cfn)

    # Create a name-getter function based on what kind of file *a* and *b* are
    get_name = gen_get_name(b, afields)

    dbed = open(BedTool._tmp(), "w")

    # keep the name, distance and feature type.
    seen_by_line = collections.defaultdict(list)
    assert len(c) > 1
    for feat in c:
        fields = feat.fields

        # Key into the dict is the reconstituted *a* feature
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

    if save_intermediates:
        dfn = 'd.bed'
    else:
        dfn = None
    d = BedTool(dbed.name).saveas(dfn)

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

    pybedtools.cleanup()

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS).failed == 0:
        main()
