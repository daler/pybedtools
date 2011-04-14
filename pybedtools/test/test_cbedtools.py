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

    def setUp(self):
        self.file = os.path.join(PATH, self.file)
        start, end, strand = 9719768, 9739768, "-"
        self.i = Interval("chr21", start, end, strand)
        self.start, self.end, self.strand = start, end, strand

    def testLengths(self):
        self.assertEqual(self.end - self.start, self.i.length)

    def testEnds(self):
        self.assertEqual(self.end, self.i.end)
        self.assertEqual(self.start, self.i.start)

    def testStrand(self):
        self.assertEqual(self.strand, self.i.strand)

    def testGetItem(self):
        """
        getitem now supports direct access to the line.
        """
        ivf = IntervalFile(self.file)
        iv = ivf.next()
        self.assert_(iv[0].startswith("chr"))
        self.assert_(iv[1].isdigit())
        self.assert_(iv[2].isdigit())

    def testGetItemSlice(self):
        """
        getitem now supports direct access to the line.
        """
        ivf = IntervalFile(self.file)
        iv = ivf.next()
        seqid, start, end = iv[0:3]
        self.assert_(start.isdigit())
        
        self.assertEqual(int(end), iv.end)
        self.assertEqual(seqid, iv.chrom)

    def testGetItemSliceNone(self):
        """
        test support for funky slices.
        """
        ivf = IntervalFile(self.file)
        iv = ivf.next()
        self.assertEqual(len(iv[:3]), 3)
        self.assertEqual(len(iv[3:3]), 0)
        self.assertEqual(len(iv[2:]), 4, iv[2:])


class IntervalFileGzTest(IntervalFileTest):
    file = "data/rmsk.hg18.chr21.small.bed.gz"


if __name__ == "__main__":
    unittest.main()


