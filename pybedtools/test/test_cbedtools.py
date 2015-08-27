#!/usr/bin/env python

import unittest
import os
from pybedtools import Interval, IntervalFile
import pybedtools
from nose.tools import assert_raises, raises
from .tfuncs import setup, teardown


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

    def testRichCmp(self):

        # be obsessive . . .
        #
        # ==
        a = Interval("chr21", 100, 200)
        b = Interval("chr21", 100, 200)
        self.assert_(a == b)
        self.assertFalse(a != b)
        self.assert_(a <= b)
        self.assert_(a >= b)
        self.assertFalse(a < b)
        self.assertFalse(a > b)

        a = Interval("chr21", 100, 100)
        b = Interval("chr21", 100, 100)
        self.assert_(a == b)
        self.assertFalse(a != b)
        self.assert_(a <= b)
        self.assert_(a >= b)
        self.assertFalse(a < b)
        self.assertFalse(a > b)


        # != because of strand
        a = Interval("chr21", 100, 200, strand='+')
        b = Interval("chr21", 100, 200, strand='-')
        self.assertFalse(a == b)
        self.assert_(a != b)
        self.assertFalse(a <= b)
        self.assertFalse(a >= b)
        self.assertFalse(a < b)
        self.assertFalse(a > b)

        # a >= b
        a = Interval("chr21", 100, 300)
        b = Interval("chr21", 100, 200)
        self.assertFalse(a == b)
        self.assert_(a != b)
        self.assertFalse(a <= b)
        self.assert_(a >= b)
        self.assertFalse(a < b)
        self.assertFalse(a > b)

        # a <= b
        a = Interval("chr21", 100, 300)
        b = Interval("chr21", 300, 300)
        self.assertFalse(a == b)
        self.assert_(a != b)
        self.assert_(a <= b)
        self.assertFalse(a >= b)
        self.assertFalse(a < b)
        self.assertFalse(a > b)


        # a <= b
        a = Interval("chr21", 100, 300)
        b = Interval("chr21", 250, 300)
        self.assertFalse(a == b)
        self.assert_(a != b)
        self.assert_(a <= b)
        self.assertFalse(a >= b)
        self.assertFalse(a < b)
        self.assertFalse(a > b)

        # a < b
        a = Interval("chr21", 100, 200)
        b = Interval("chr21", 201, 300)
        self.assertFalse(a == b)
        self.assert_(a != b)
        self.assert_(a <= b)
        self.assertFalse(a >= b)
        self.assert_(a < b)
        self.assertFalse(a > b)

        # a > b
        a = Interval("chr21", 201, 300)
        b = Interval("chr21", 100, 200)
        self.assertFalse(a == b)
        self.assert_(a != b)
        self.assertFalse(a <= b)
        self.assert_(a >= b)
        self.assertFalse(a < b)
        self.assert_(a > b)

        # a != b
        a = Interval("none", 1, 100)
        b = Interval("chr21", 1, 100)
        self.assertFalse(a == b)
        self.assert_(a != b)
        self.assertFalse(a <= b)
        self.assertFalse(a >= b)
        self.assertFalse(a < b)
        self.assertFalse(a > b)

        # nested should raise NotImplementedError
        a = Interval("chr21", 100, 200)
        b = Interval("chr21", 50, 300)
        self.assertRaises(NotImplementedError, a.__eq__, b)
        self.assertRaises(NotImplementedError, a.__ne__, b)
        self.assertRaises(NotImplementedError, a.__le__, b)
        self.assertRaises(NotImplementedError, a.__ge__, b)
        self.assertRaises(NotImplementedError, a.__lt__, b)
        self.assertRaises(NotImplementedError, a.__gt__, b)


class IntervalTest(unittest.TestCase):
    file = "data/rmsk.hg18.chr21.small.bed.gz"
    chrpos = 0
    startpos = 1
    stoppos = 2
    fieldcount = 6

    def setUp(self):
        self.file = os.path.join(PATH, self.file)
        start, end, strand = 9719768, 9739768, "-"
        self.i = Interval("chr21", start, end, strand=strand)
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
        iv = next(ivf)
        self.assert_(iv[self.chrpos].startswith("chr"))
        self.assert_(iv[self.startpos].isdigit())
        self.assert_(iv[self.startpos].isdigit())

    def testGetItemNegative(self):
        "test negative indexes to feature."
        ivf = IntervalFile(self.file)
        iv = next(ivf)
        self.assert_(iv[-self.fieldcount+self.chrpos].startswith("chr"), iv[-self.fieldcount+self.chrpos])
        self.assert_(iv[-self.fieldcount+self.startpos].isdigit(), iv[-self.fieldcount+self.startpos])
        self.assert_(iv[-self.fieldcount+self.stoppos].isdigit())

    def testGetItemSlice(self):
        "getitem now supports direct access to the line."
        ivf = IntervalFile(self.file)
        iv = next(ivf)
        seqid, = iv[self.chrpos:self.chrpos+1]
        start, end = iv[self.startpos:self.stoppos+1]
        self.assert_(start.isdigit())

        self.assertEqual(int(end), iv.end)
        self.assertEqual(seqid, iv.chrom)

    def testGetItemSliceNone(self):
        " test support for funky slices."
        ivf = IntervalFile(self.file)
        iv = next(ivf)
        self.assertEqual(len(iv[:3]), 3)
        self.assertEqual(len(iv[3:3]), 0)
        self.assertEqual(len(iv[2:]), self.fieldcount-2, iv[2:])
        
        print(len(iv.fields), iv.fields)
        self.assertRaises(IndexError, lambda x: iv[x], self.fieldcount+1)

    def testGetItemString(self):
        ivf = IntervalFile(self.file)
        iv = next(ivf)
        self.assertEqual(iv['chrom'], iv.chrom)
        self.assertEqual(iv['start'], iv.start)
        self.assertEqual(iv['end'], iv.end)

    def testSetItemString(self):
        ivf = IntervalFile(self.file)
        iv = next(ivf)
        iv['chrom'] = 'fake'
        self.assertEqual(iv['chrom'], 'fake')
        self.assertEqual(iv.chrom, 'fake')

    #TODO: need some work on getting and setting before running these
    def testSetItem(self):
        ivf = IntervalFile(self.file)
        iv = next(ivf)
        iv.chrom = 'chrfake'
        print(iv.fields)
        self.assertEqual(iv['chrom'], 'chrfake')
        self.assertEqual(iv.chrom, 'chrfake')

    def testSetAttrs(self):
        ivf = IntervalFile(self.file)
        iv = next(ivf)
        if iv.file_type != 'gff':
            iv.attrs['a'] = 'b'
            self.assertRaises(ValueError, str, iv)
            return
        iv.attrs['ID'] = 'fake'
        iv.attrs['field0'] = 'asdf'
        self.assertEqual(str(iv.attrs), iv[8])
        self.assert_('field0=asdf' in iv[8])
        self.assert_('ID=fake' in iv[8])

    def testAppend(self):
        ivf = IntervalFile(self.file)
        iv = next(ivf)
        print(iv.fields)
        iv.append('asdf')
        print(iv)
        self.assertEqual(iv[-1], 'asdf')

    def testName(self):
        ivf = IntervalFile(self.file)
        iv = next(ivf)
        iv.name = "bart simpson"
        self.assertEqual(iv.name, "bart simpson")
        if iv.file_type == "gff":
            self.assert_("bart" in iv.fields[8])

    def testStart(self):
        ivf = IntervalFile(self.file)
        iv = next(ivf)
        orig_string = str(iv)
        orig_start = iv.start
        iv.start = orig_start
        second_string = str(iv)
        second_start = iv.start
        iv.start = second_start
        print('   orig:', '(start=%s)'%orig_start, orig_string)
        print(' second:', '(start=%s)'%second_start, second_string)
        print('current:', '(start=%s)'%iv.start, str(iv))
        self.assert_(orig_start == second_start == iv.start)
        self.assert_(orig_string == second_string == str(iv))


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
        self.i = Interval("chr1", start, end, strand=strand)
        self.start, self.end, self.strand = start, end, strand

    # Overwrite IntervalTest.testStart
    def testStart(self):
        ivf = IntervalFile(self.file)
        iv = next(ivf)
        orig_string = str(iv)

        # 0-based.
        orig_start = iv.start

        # Setting .start always sets 0-based coord.
        iv.start = orig_start

        # But for GFF setting .start should also make the .fields[3] the GFF
        # 1-based coord
        assert iv.start == int(iv.fields[3])-1

        second_string = str(iv)
        second_start = iv.start
        iv.start = second_start

        # Check .start and .fields[3] internal consistency again
        assert iv.start == int(iv.fields[3])-1

        print('   orig:', '(start=%s)'%orig_start, orig_string)
        print(' second:', '(start=%s)'%second_start, second_string)
        print('current:', '(start=%s)'%iv.start, str(iv))
        self.assert_(orig_start == second_start == iv.start)
        self.assert_(orig_string == second_string == str(iv))


def test_zero_length_regression():
    i = Interval(chrom='chrA', start=1, end=10, name='.', score='0', strand='+')
    iminus = Interval(chrom='chrA', start=1, end=10, name='.', score='0', strand='-')
    m = pybedtools.BedTool('chrA 2 2 . 0 +', from_string=True).intervals

    assert len(m.all_hits(i)) == 1, m.all_hits(i)
    assert len(m.all_hits(iminus)) == 1
    assert len(m.all_hits(i, same_strand=False)) == 1
    assert len(m.all_hits(iminus, same_strand=False)) == 1
    assert len(m.all_hits(i, same_strand=True)) == 1
    assert len(m.all_hits(iminus, same_strand=True)) == 0

    assert m.any_hits(i) == 1
    assert m.any_hits(iminus) == 1
    assert m.any_hits(i, same_strand=True) == 1
    assert m.any_hits(i, same_strand=False) == 1
    assert m.any_hits(iminus, same_strand=True) == 0
    assert m.any_hits(iminus, same_strand=False) == 1

    i2 = Interval(chrom='chrA', start=999, end=9999)
    assert len(m.all_hits(i2)) == 0
    assert len(m.all_hits(i2, same_strand=True)) == 0

    a = pybedtools.BedTool('chr1 3 3', from_string=True)
    print([len(i) for i in a.intervals])
    result = a.intervals.all_hits(Interval('chr1', 3, 3))
    assert len(result) == 1, result

def test_issue_123():
    bed = """\
chr1    1       1       Region_A        2       +
chr1    2       2       Region_B        1       +
chr1    3       3       Region_C        3       +
chr1    4       5
""".strip()

    b = pybedtools.BedTool(bed, from_string=True).intervals

    result = b.all_hits(Interval('chr1', 3, 3))
    assert len(result) == 1, len(result)

def test_missing_files():
    """
    previously this would crash the interpreter due to an exit(1) call in
    bedFile.cpp
    """
    a = pybedtools.BedTool('chrA 1 10', from_string=True).saveas("this_file_should_raise_BEDTools_Error")
    result = list(iter(a))
    os.unlink(a.fn)
    from pybedtools.cbedtools import BedToolsFileError
    def crashes():
        list(iter(a))

    assert_raises(BedToolsFileError, crashes)

if __name__ == "__main__":
    unittest.main()
    pybedtools.cleanup(remove_all=True)


