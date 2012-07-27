from collections import defaultdict
import matplotlib
from matplotlib import collections
from matplotlib import pyplot as plt
import numpy as np
import pybedtools


class Track(collections.PolyCollection):
    def __init__(self, features, chrom=None, ybase=0, yheight=1,
            visibility='dense', stranded=True, **kwargs):
        """
        Subclass of matplotlib's PolyCollection that can be added to an Axes.

        :param features:
            Can be an existing BedTool, or anything than can be used to create
            a BedTool (e.g., a filename or a generator of Interval objects)

        :param ybase:
            y-coord of the bottom edge of the track (in data coordinates)

        :param yheight:
            How high each feature will be, in data coordinates

        :param visibility:
            Mimics the settings in the UCSC Genome Browser:

            * "dense" is the default; overlapping features can be seen if you
              set alpha < 1.

            * "squish" prevents adjacent features from overlapping.  This keeps
              `yheight` for all features, so if you have a lot of features
              piling up, the track will be a lot higher on the y-axis than
              `yheight`.

        :param stranded:
            If boolean and True, will draw arrrow-shaped features to indicate
            direction (where the point is 10% of the total gene length)

            If a dictionary, map strands to colors, e.g., {'+': 'r', '-': 'b'}.

        :param kwargs:
            Additional keyword args are passed to
            matplotlib.collections.PolyCollection.

        Notes:

        After creating a track, use the `ymax` attribute to get the max y-value
        used in the track -- useful if you've created a "squish" track but
        would like to stack another track on top, and need to calculate what
        the new Track's `ybase` should be.

        The returned PolyCollection will have the `features` attribute, which
        contains the BedTool it was created from -- so you can write callback
        functions for event handling, e.g.::

            def callback(event):
                '''
                prints the feature's line when clicked in the plot
                '''
                coll = event.artist
                for i in event.ind:
                    print coll.features[i]

            fig.canvas.mpl_connect('on_pick', callback)


        >>> a = pybedtools.example_bedtool('a.bed')
        >>> track = Track(a, alpha=0.5, picker=5)
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111)
        >>> ax.add_collection(track) #doctest: +ELLIPSIS
        <pybedtools.contrib.plotting.Track object at 0x...>
        >>> limits = ax.axis('tight')
        """
        if isinstance(features, pybedtools.BedTool)\
                and isinstance(features.fn, basestring):
            self.features = features
        else:
            self.features = pybedtools.BedTool(features).saveas()
        self._visibility = visibility
        self._ybase = ybase
        self._yheight = yheight
        self.stranded = stranded
        self._check_stranded_dict()
        facecolors = self._colors()
        kwargs.update(dict(facecolors=facecolors))
        collections.PolyCollection.__init__(
                self, verts=self._get_verts(), **kwargs)

    def _shape(self, feature, ybase, yheight):
        if self.stranded and not isinstance(self.stranded, dict):
            offset = len(feature) * 0.1
            if feature.strand == '-':
                return [
                        (feature.stop, ybase),
                        (feature.stop, ybase + yheight),
                        (feature.start + offset, ybase + yheight),
                        (feature.start, ybase + yheight * 0.5),
                        (feature.start + offset, ybase)
                        ]

            elif feature.strand == '+':
                return [
                        (feature.start, ybase),
                        (feature.start, ybase + yheight),
                        (feature.stop - offset, ybase + yheight),
                        (feature.stop, ybase + yheight * 0.5),
                        (feature.stop - offset, ybase)
                        ]
        return [
                (feature.start, ybase),
                (feature.start, ybase + yheight),
                (feature.stop, ybase + yheight),
                (feature.stop, ybase)
                ]

    def _get_verts(self):
        verts = []

        if self._visibility == 'dense':
            for feature in self.features:
                verts.append(self._shape(feature, self._ybase, self._yheight))
            self.ymax = self._ybase + self._yheight

        if self._visibility == 'squish':
            # Using "squish" mode will create multiple "strata" of features.
            # The stack keeps track of the end coord of the longest feature in
            # each strata
            #
            # Reasonably efficient -- <2s to plot 15K multiply-overlapping
            # features
            stack = []
            ybase = self._ybase
            self.ymax = self._ybase + self._yheight
            for feature in self.features:
                ybase = None
                for i, s in enumerate(stack):
                    if feature.start > s:
                        ybase = self._ybase + i * self._yheight
                        stack[i] = feature.stop
                        break
                if ybase is None:
                    ybase = self._ybase + len(stack) * self._yheight
                    stack.append(feature.end)
                verts.append(self._shape(feature, ybase, self._yheight))
            self.ymax = self._ybase + len(stack) * self._yheight

        return verts

    def _check_stranded_dict(self):
        if not isinstance(self.stranded, dict):
            return True
        if '+' not in self.stranded:
            raise ValueError('stranded dict "%s" does not have required '
                    'key "+"' % self.stranded)
        if '-' not in self.stranded:
            raise ValueError('stranded dict "%s" does not have required '
                    'key "-"' % self.stranded)
        return True

    def _colors(self):
        if not isinstance(self.stranded, dict):
            return None
        colors = []
        for feature in self.features:
            try:
                colors.append(self.stranded[feature.strand])
            except KeyError:
                raise KeyError('strand color dict "%s" does not have a key '
                        'for strand "%s"' % (self.stranded, feature.strand))
        return colors

    def get_xlims(self, ax):
        """
        Needs `ax` to convert to transData coords
        """
        bb = self.get_datalim(ax.transData)
        return (bb.xmin, bb.xmax)

    @property
    def midpoint(self):
        return self._ybase + (self.ymax - self._ybase) / 2.0


def binary_heatmap(bts, names, plot=True):
    """
    Plots a "binary heatmap", showing the results of a multi-intersection.

    Each row is a different genomic region found in at least one of the input
    BedTools; each column represents a different file.  Black indicates whether
    a feature was found at that particular site.  Rows with black all the way
    across indicates that all features were colocalized at those sites.

    `bts` is an iterable of BedTool objects or filenames; `names` is a list of
    labels to use in the plot and is exactly the same length as `bts`.

    If `plot=True`, then plot the sorted, labeled matrix with matplotlib.

    Returns (summary, m) where `summary` is a dictionary summarizing the
    results and `m` is the sorted NumPy array.  See source for further details.
    """
    # Be flexible about input types
    _bts = []
    for bt in bts:
        if isinstance(bt, pybedtools.BedTool):
            if not isinstance(bt.fn, basestring):
                bt = bt.saveas()
            _bts.append(bt.fn)
        elif isinstance(bt, basestring):
            _bts.append(bt)

    # Do the multi-intersection.
    results = pybedtools.BedTool().multi_intersect(
            i=_bts,
            names=names,
            cluster=True)

    # If 4 files were provided with labels 'a', 'b', 'c', and 'd, each line
    # would look something like:
    #
    #   chr2L    65716    65765    4    a,b,c,d    1    1    1    1
    #   chr2L    71986    72326    1    c          0    0    1    0
    #
    # The last four columns will become the matrix; save the class labels (5th
    # column) for a printed out report
    d = defaultdict(int)
    m = []
    for item in results:
        cls = item[4]
        d[cls] += 1
        m.append(item[5:])

    m = np.array(m, dtype=int)
    ind = sort_binary_matrix(m)

    if plot:
        # Plot and label it
        fig = plt.figure(figsize=(3, 10))
        ax = fig.add_subplot(111)

        # matplotlib.cm.binary: 1 = black, 0 = white; force origin='upper' so
        # that array's [0,0] is in the upper left corner.
        mappable = ax.imshow(m[ind], aspect='auto', interpolation='nearest',
                cmap=matplotlib.cm.binary, origin='upper')
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=90)
        fig.subplots_adjust(left=0.25)

    return d, m


def sort_binary_matrix(m):
    """
    Performs a column-weighted sort on a binary matrix, returning the new index
    """
    # To impart some order in the matrix, give columns increasingly higher
    # weights...
    weights = [2 ** i for i in range(1, m.shape[1] + 1)[::-1]]

    # ...then create scores...
    score_mat = m * weights

    # ...and re-sort the matrix based on row sums (reversed so that highest
    # scores are on top)
    ind = np.argsort(score_mat.sum(axis=1))[::-1]
    return ind


def binary_summary(d):
    """
    Convenience function useful printing the results from binary_heatmap().
    """
    s = []
    for item in sorted(d.items(), key=lambda x: x[1], reverse=True):
        s.append('%s : %s' % (item))
    return '\n'.join(s)

if __name__ == "__main__":
    bts = [
            pybedtools.example_bedtool('BEAF_Kc_Bushey_2009.bed'),
            pybedtools.example_bedtool('CTCF_Kc_Bushey_2009.bed'),
            pybedtools.example_bedtool('Cp190_Kc_Bushey_2009.bed'),
            pybedtools.example_bedtool('SuHw_Kc_Bushey_2009.bed'),
        ]
    names = ['BEAF', 'CTCF', 'Cp190', 'Su(Hw)']

    #bts = [
    #        pybedtools.example_bedtool('a.bed'),
    #        pybedtools.example_bedtool('b.bed')]
    #names = ['a','b']
    d, m = binary_heatmap(bts, names)
    print binary_summary(d)
    plt.show()
