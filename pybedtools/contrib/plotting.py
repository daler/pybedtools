import matplotlib
from matplotlib import collections
from matplotlib import pyplot as plt
import pybedtools


class Track(collections.PolyCollection):
    def __init__(self, features, chrom=None, ybase=0, yheight=1,
            visibility='dense', **kwargs):
        """
        Subclass of matplotlib's PolyCollection that can be added to an Axes.

        `features` can be an existing BedTool, or anything than can be used to
        create a BedTool (e.g., a filename or a generator of Interval objects)

        `ybase` is the y-coord of the bottom edge of the track.

        `yheight` is how high each feature will be.

        `visibility` mimics the settings in the UCSC Genome Browser:

            * "dense" is the default; overlapping features can be seen if you
              set alpha < 1.

            * "squish" prevents adjacent features from overlapping.  This keeps
              `yheight` for all features, so if you have a lot of features
              piling up, the track will be a lot higher on the y-axis than
              `yheight`.


        `**kwargs` are passed to matplotlib.collections.PolyCollection.

        Use Track.ymax to get the max y-value used in the track -- useful if
        you've created a "squish" track but would like to stack another track
        on top, and need to calculate what the new Track's `ybase` should be.

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
        >>> ax.add_collection(track)
        >>> ax.axis('tight')
        """
        if isinstance(features, pybedtools.BedTool)\
                and isinstance(features.fn, basestring):
            self.features = features
        else:
            self.features = pybedtools.BedTool(features).saveas()
        self._visibility = visibility
        self._ybase = ybase
        self._yheight = yheight
        collections.PolyCollection.__init__(
                self, verts=self._get_verts(), **kwargs)

    def _get_verts(self):
        verts = []

        if self._visibility == 'dense':
            for feature in self.features:
                verts.append([
                    (feature.start, self._ybase),
                    (feature.start, self._ybase + self._yheight),
                    (feature.stop, self._ybase + self._yheight),
                    (feature.stop, self._ybase),
                    ])
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
                verts.append([
                    (feature.start, ybase),
                    (feature.start, ybase + self._yheight),
                    (feature.stop, ybase + self._yheight),
                    (feature.stop, ybase),
                    ])
            self.ymax = self._ybase + len(stack) * self._yheight

        return verts

    @property
    def midpoint(self):
        return self._ybase + (self.ymax - self._ybase) / 2.0
