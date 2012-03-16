import matplotlib
from matplotlib import collections
from matplotlib import pyplot as plt
import pybedtools


class Track(object):
    def __init__(self, bed):
        """
        `bed` can be an existing BedTool, or anything than can be used to
        create a BedTool (e.g., a filename or a generator of Interval objects)
        >>> a = pybedtools.example_bedtool('a.bed')
        >>> track = Track(a)

        >>> a_fn = pybedtools.example_filename('a.bed')
        >>> track = Track(a_fn)
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111)
        >>> ax.add_collection(track.collection(alpha=0.5, picker=5))
        >>> ax.axis('tight')
        """
        if isinstance(bed, pybedtools.BedTool)\
                and isinstance(bed.fn, basestring):
            self.bed = bed
        else:
            self.bed = pybedtools.BedTool(bed).saveas()

    def collection(self, ybase=0, yheight=1,
            **kwargs):
        """
        `ybase` is the y-coord of the bottom edge of each feature.

        `yheight` is how high each feature will be.

        `**kwargs` are passed to matplotlib.collections.PolyCollection.

        The returned PolyCollection will have the `features` attribute, which
        contains the BedTool it was created from -- so you can write callback
        functions for event handling, e.g.::

            def callback(event):
                coll = event.artist
                for i in event.ind:
                    print coll.features[i]

            fig.canvas.mpl_connect('on_pick', callback)

        """
        verts = []
        for feature in self.bed:
            verts.append([
                (feature.start, ybase),
                (feature.start, ybase + yheight),
                (feature.stop, ybase + yheight),
                (feature.stop, ybase),
                ])
        coll = matplotlib.collections.PolyCollection(verts, **kwargs)
        coll.features = self.bed
        return coll
