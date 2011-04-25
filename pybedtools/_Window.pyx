# cython: profile=True

import os
from collections import deque

cdef class Window(object):
    cdef public object iterable
    cdef public object low_reads
    cdef public object high_reads
    cdef public int windowsize
    cdef int center
    cdef int left_edge
    cdef int right_edge
    cdef str chrom
    cdef public int debug
    cdef public object buffered_read
    cdef int START

    def __init__(self, iterable, windowsize=100, debug=0):
        """
        Constructor:

            Window(iterable, windowsize=100, debug=0)

        Moving window over an *iterable* of features (e.g., IntervalFile(fn)) of
        size *windowsize*.  Use *debug=1* to see all sorts of output for
        double-checking.

        The resulting Window instance can be iterated over.  Each iteration
        returns a tuple of::

            (center, low_reads, high_reads)

        where *center* is the current center of the window; *low_reads* is a
        deque of reads that includes the center and everything lower than it
        that fits within the window; and *high_reads* is a deque of reads that
        includes everything higher within the window.

        The strategy is to hold one read as the "centered read", which is
        currently in focus.  Reads are checked to see if they fit within the
        window centered on this read.  There is always a buffered read, which
        is last read taken from the iterable.  If the buffered read doesn't fit
        in the window, it remains the buffered read until the current window is
        returned.  It will continue to remain the buffered read (and no more
        reads will be taken from the iterable) until it fits within the current
        window.

        The window is implemented in two parts, a low_reads and a high_reads
        part.

        The next window's center jumps to the next available read position,
        rather than the next available bp.  This can greatly save on running
        time.  The next available read position will typically be the first
        item in the high_reads deque.

        """
        self.iterable = iterable
        self.windowsize = windowsize
        self.left_edge = 0
        self.right_edge = 0
        self.debug = debug

        # Here we pull the first thing from the iterable to set up the various
        # attributes
        first_read = self.iterable.next()
        self.chrom = first_read.chrom
        first_start_pos = first_read.start
        self.left_edge = first_start_pos - self.windowsize/2
        self.right_edge = self.left_edge + self.windowsize
        self.center = first_start_pos
        self.buffered_read = first_read
        self.high_reads = deque()
        self.low_reads = deque([self.buffered_read])
        self.START = 1

    cdef int accumulate_reads(self) except -1:
        """
        Fill up the window surrounding the currently-centered read.
        """
        if self.debug:
            print 'appending:\n\t',

        while True:


            # Need to short-circuit if starting, cause we've already filled
            # buffered_read
            if self.START:
                self.START = 0
                self.buffered_read = self.iterable.next()
                continue

            if self.buffered_read.chrom != self.chrom:
                if self.debug:
                    print 'new chrom -- %s' % self.buffered_read.chrom
                break

            # While accumulating, the only time low_reads will fill up is if
            # they are duplicates of the currently-centered read
            if self.buffered_read.start == self.center:

                if self.debug:
                    print self.buffered_read.start,

                self.low_reads.append(self.buffered_read)

            # Otherwise, if it's within the window then it's added to
            # high_reads.
            elif self.buffered_read.start < self.right_edge:

                if self.debug:
                    print  self.buffered_read.start,

                self.high_reads.append(self.buffered_read)

            else:
                break

            # The positioning of this is important -- we only get a new
            # buffered read if the last buffered read has been treated --
            # either added to low_reads or high_reads
            self.buffered_read = self.iterable.next()

        if self.debug:
            print

        return 0

    cdef int trim(self):
        """
        Trims reads off window edges, which is basically just shifting the
        window.
        """

        # If there is nothing in the high reads, then use the current buffered
        # read as the center.
        if len(self.high_reads) == 0:
            self.center = self.buffered_read.start
            self.chrom = self.buffered_read.chrom
            self.left_edge = self.center - self.windowsize/2
            self.right_edge = self.center + self.windowsize/2

        # Otherwise, use the next read in the high_reads deque
        else:
            self.chrom = self.high_reads[0].chrom
            self.center = self.high_reads[0].start
            self.left_edge = self.center - self.windowsize/2
            self.right_edge = self.center + self.windowsize/2

        # Now that the center point has been reset, remove reads from low_reads
        # list that no longer fit in the window
        if self.debug:
            print 'removed:', 
        while True:

            # Must be a better way to do this other than popping it off and
            # then back on again if it's in range, though the appendleft will
            # only happen during one (i.e. the last) time through the loop
            try:
                popped = self.low_reads.popleft()
                if (popped.start < self.left_edge) or (popped.chrom != self.buffered_read.chrom):
                    if self.debug:
                        print popped.start,
                    continue
                else:
                    self.low_reads.appendleft(popped)
                    break

            # If there's nothing left in the low_reads, then stop removing
            except IndexError:
                break

        # Next we remove any additional reads that are duplicates of the
        # centered read and add these to low_reads
        while True: 
            try:
                popped = self.high_reads.popleft()
                if popped.start == self.center:
                    self.low_reads.append(popped)
                else:
                    self.high_reads.appendleft(popped)
                    break
            except IndexError:
                break

        # Run accumulator again to see if we can add the current buffered read
        # and/or any additional reads to the window.
        #self.accumulate_reads()
        if self.debug:
            print 

        return 0

    def __iter__(self):
        return self

    def __next__(self):

        if not self.START:
            # This moves the window...
            self.trim()

        # First we accumulate reads
        self.accumulate_reads()

        if self.debug:
            print 'chrom        :', self.chrom
            print 'left         :', self.left_edge
            print 'center       :', self.center
            print 'right        :', self.right_edge
            print 'low contents :', [i.start for i in self.low_reads]
            print 'high contents:', [i.start for i in self.high_reads]
            print 'buffer       :', self.buffered_read.start

        return self.center, self.low_reads, self.high_reads


