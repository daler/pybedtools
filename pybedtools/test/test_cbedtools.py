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

    def setUp(self):
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



class IntervalFileGzTest(IntervalFileTest):
      file = "data/rmsk.hg18.chr21.small.bed.gz"


if __name__ == "__main__":
    unittest.main()


