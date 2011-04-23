#!/usr/bin/env python

import unittest
import os
from pybedtools import Interval, IntervalFile
PATH = os.path.dirname(__file__)

class IntervalFileTest(unittest.TestCase):
    file = "data/rmsk.hg18.chr21.small.bed"
    def setUp(self):
        self.file = os.path.join(PATH, self.file)
        self.bed = IntervalFile(self.file)

    def testFileType(self):
        self.assert_(self.bed.file_type == "bed", (self.bed.file_type, self.file))

        gff = os.path.join(PATH, "data/c.gff")
        i = IntervalFile(gff)
        self.assert_(i.file_type == "gff", (i.file_type, gff))

    def testOverlaps(self):
        i    = Interval("chr21", 9719768, 9739768)
        hits = self.bed.all_hits(i)
        self.assertEqual(len(hits), 8)
        for hit in hits:
            self.assert_(hit.start <= 9739768 and hit.end >= 9719768)

    def testStrands(self):
        i = Interval("chr21", 9719768, 9739768, "+")
        hits = self.bed.all_hits(i, same_strand=True)
        for hit in hits:
            self.assert_(hit.strand == '+')

        i = Interval("chr21", 9719768, 9739768, "-")
        hits = self.bed.all_hits(i, same_strand=True)
        for hit in hits:
            self.assert_(hit.strand == '-')

class IntervalTest(unittest.TestCase):
    file = "data/rmsk.hg18.chr21.small.bed.gz"
    chrpos = 0
    startpos = 1
    stoppos = 2
    fieldcount = 6

    def setUp(self):
        self.file = os.path.join(PATH, self.file)
        start, end, strand = 9719768, 9739768, "-"
        self.i = Interval("chr21", start, end, strand)
        self.start, self.end, self.strand = start, end, strand

    def testLengths(self):
        self.assertEqual(self.end - self.start, self.i.length)
        self.assertEqual(len(self.i), self.i.length)

    def testEnds(self):
        self.assertEqual(self.end, self.i.end)
        self.assertEqual(self.start, self.i.start)

    def testStrand(self):
        self.assertEqual(self.strand, self.i.strand)

    def testGetItem(self):
        "getitem now supports direct access to the line."
        ivf = IntervalFile(self.file)
        iv = ivf.next()
        self.assert_(iv[self.chrpos].startswith("chr"))
        self.assert_(iv[self.startpos].isdigit())
        self.assert_(iv[self.startpos].isdigit())

    def testGetItemNegative(self):
        "test negative indexes to feature."
        ivf = IntervalFile(self.file)
        iv = ivf.next()
        self.assert_(iv[-self.fieldcount+self.chrpos].startswith("chr"), iv[-self.fieldcount+self.chrpos])
        self.assert_(iv[-self.fieldcount+self.startpos].isdigit(), iv[-self.fieldcount+self.startpos])
        self.assert_(iv[-self.fieldcount+self.stoppos].isdigit())

    def testGetItemSlice(self):
        "getitem now supports direct access to the line."
        ivf = IntervalFile(self.file)
        iv = ivf.next()
        seqid, = iv[self.chrpos:self.chrpos+1]
        start, end = iv[self.startpos:self.stoppos+1]
        self.assert_(start.isdigit())

        self.assertEqual(int(end), iv.end)
        self.assertEqual(seqid, iv.chrom)

    def testGetItemSliceNone(self):
        " test support for funky slices."
        ivf = IntervalFile(self.file)
        iv = ivf.next()
        self.assertEqual(len(iv[:3]), 3)
        self.assertEqual(len(iv[3:3]), 0)
        self.assertEqual(len(iv[2:]), self.fieldcount-2, iv[2:])
        
        print len(iv.fields), iv.fields
        self.assertRaises(IndexError, lambda x: iv[x], self.fieldcount+1)

    def testGetItemString(self):
        ivf = IntervalFile(self.file)
        iv = ivf.next()
        self.assertEqual(iv['chrom'], iv.chrom)
        self.assertEqual(iv['start'], iv.start)
        self.assertEqual(iv['end'], iv.end)

    def testSetItemString(self):
        ivf = IntervalFile(self.file)
        iv = ivf.next()
        iv['chrom'] = 'fake'
        self.assertEqual(iv['chrom'], 'fake')
        self.assertEqual(iv.chrom, 'fake')

    #TODO: need some work on getting and setting before running these
    def testSetItem(self):
        ivf = IntervalFile(self.file)
        iv = ivf.next()
        iv.chrom = 'chrfake'
        print iv.fields
        self.assertEqual(iv['chrom'], 'chrfake')
        self.assertEqual(iv.chrom, 'chrfake')

    def testAppend(self):
        ivf = IntervalFile(self.file)
        iv = ivf.next()
        print iv.fields
        iv.append('asdf')
        print iv
        self.assertEqual(iv[-1], 'asdf')

    def testName(self):
        ivf = IntervalFile(self.file)
        iv = ivf.next()
        iv.name = "bart simpson"
        self.assertEqual(iv.name, "bart simpson")
        if iv.file_type == "gff":
            self.assert_("bart" in iv.fields[8])

class IntervalFileGzTest(IntervalFileTest):
    file = "data/rmsk.hg18.chr21.small.bed.gz"

class IntervalFileGFFTest(IntervalTest):
    file = 'data/d.gff'
    chrpos = 0
    startpos = 3
    stoppos = 4
    fieldcount = 9

    def setUp(self):
        self.file = os.path.join(PATH, self.file)
        start, end, strand = 1, 100, "+"
        self.i = Interval("chr1", start, end, strand)
        self.start, self.end, self.strand = start, end, strand



if __name__ == "__main__":
    unittest.main()


