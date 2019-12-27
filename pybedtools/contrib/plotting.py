import os
from collections import defaultdict
import matplotlib
from matplotlib import collections
from matplotlib import pyplot as plt
import numpy as np
import pybedtools
import six


class Track(collections.PolyCollection):
    def __init__(
        self,
        features,
        chrom=None,
        ybase=0,
        yheight=1,
        visibility="dense",
        stranded=True,
        **kwargs
    ):
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
        >>> track = pybedtools.contrib.plotting.Track(a, alpha=0.5, picker=5)
        >>> import matplotlib.pyplot as plt
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111)
        >>> ax.add_collection(track) #doctest: +ELLIPSIS
        <pybedtools.contrib.plotting.Track object at 0x...>
        >>> limits = ax.axis('tight')
        """
        if isinstance(features, pybedtools.BedTool) and isinstance(
            features.fn, six.string_types
        ):
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
        collections.PolyCollection.__init__(self, verts=self._get_verts(), **kwargs)

    def _shape(self, feature, ybase, yheight):
        if self.stranded and not isinstance(self.stranded, dict):
            offset = len(feature) * 0.1
            if feature.strand == "-":
                return [
                    (feature.stop, ybase),
                    (feature.stop, ybase + yheight),
                    (feature.start + offset, ybase + yheight),
                    (feature.start, ybase + yheight * 0.5),
                    (feature.start + offset, ybase),
                ]

            elif feature.strand == "+":
                return [
                    (feature.start, ybase),
                    (feature.start, ybase + yheight),
                    (feature.stop - offset, ybase + yheight),
                    (feature.stop, ybase + yheight * 0.5),
                    (feature.stop - offset, ybase),
                ]
        return [
            (feature.start, ybase),
            (feature.start, ybase + yheight),
            (feature.stop, ybase + yheight),
            (feature.stop, ybase),
        ]

    def _get_verts(self):
        verts = []

        if self._visibility == "dense":
            for feature in self.features:
                verts.append(self._shape(feature, self._ybase, self._yheight))
            self.ymax = self._ybase + self._yheight

        if self._visibility == "squish":
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
        if "+" not in self.stranded:
            raise ValueError(
                'stranded dict "%s" does not have required ' 'key "+"' % self.stranded
            )
        if "-" not in self.stranded:
            raise ValueError(
                'stranded dict "%s" does not have required ' 'key "-"' % self.stranded
            )
        return True

    def _colors(self):
        if not isinstance(self.stranded, dict):
            return None
        colors = []
        for feature in self.features:
            try:
                colors.append(self.stranded[feature.strand])
            except KeyError:
                raise KeyError(
                    'strand color dict "%s" does not have a key '
                    'for strand "%s"' % (self.stranded, feature.strand)
                )
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


class BinaryHeatmap(object):
    """
    Class-based version of the `binary_heatmap` function for more flexibility.
    """

    def __init__(self, bts, names):
        self.bts = bts
        self.names = names

        # Be flexible about input types
        _bts = []
        for bt in bts:
            if isinstance(bt, pybedtools.BedTool):
                if not isinstance(bt.fn, six.string_types):
                    bt = bt.saveas()
                _bts.append(bt.fn)
            elif isinstance(bt, six.string_types):
                _bts.append(bt)

        # Do the multi-intersection.
        self.results = pybedtools.BedTool().multi_intersect(
            i=_bts, names=names, cluster=True
        )

        # If 4 files were provided with labels 'a', 'b', 'c', and 'd, each line
        # would look something like:
        #
        #   chr2L    65716    65765    4    a,b,c,d    1    1    1    1
        #   chr2L    71986    72326    1    c          0    0    1    0
        #
        # The last four columns will become the matrix; save the class labels (5th
        # column) for a printed out report
        self.class_counts = defaultdict(int)
        _classified_intervals = defaultdict(list)
        self.matrix = []
        for item in self.results:
            cls = item[4]
            self.class_counts[cls] += 1
            self.matrix.append(item[5:])
            _classified_intervals[cls].append(item)

        self.classified_intervals = {}
        for k, v in list(_classified_intervals.items()):
            self.classified_intervals[k] = pybedtools.BedTool(v)

        self.matrix = np.array(self.matrix, dtype=int)
        self.sort_ind = sort_binary_matrix(self.matrix)

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure(figsize=(3, 10))
            ax = fig.add_subplot(111)
        # matplotlib.cm.binary: 1 = black, 0 = white; force origin='upper' so
        # that array's [0,0] is in the upper left corner.
        mappable = ax.imshow(
            self.matrix[self.sort_ind],
            aspect="auto",
            interpolation="nearest",
            cmap=matplotlib.cm.binary,
            origin="upper",
        )
        ax.set_xticks(list(range(len(self.names))))
        ax.set_xticklabels(self.names, rotation=90)
        if ax is None:
            fig.subplots_adjust(left=0.25)
        return ax


def binary_heatmap(bts, names, plot=True, cluster=True):
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
    bh = BinaryHeatmap(bts=bts, names=names)
    if plot:
        bh.plot()

    return bh.class_counts, bh.matrix


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
    for item in sorted(list(d.items()), key=lambda x: x[1], reverse=True):
        s.append("%s : %s" % (item))
    return "\n".join(s)


class TrackCollection(object):
    def __init__(self, config, yheight=1, figsize=None, padding=0.1):
        """
        Handles multiple tracks on the same figure.

        :param config:
            A list of tuples that configures tracks.

            Each tuple contains a filename, BedTool object, or other
            iterable of pybedtools.Interval objects and a dictionary of
            keyword args that will be used to create a corresponding Track
            object, e.g.::

                [
                    ('a.bed',
                        dict(color='r', alpha=0.5, label='a')),
                    (BedTool('a.bed').intersect('b.bed'),
                        dict(color='g', label='b')),
                ]

            In this dictionary, do not specify `ybase`, since that will be
            handled for you.  Also do not specify `yheight` in these
            dictionaries -- `yheight` should be provided as a separate kwarg to
            so that the `padding` kwarg works correctly.

        :param figsize:
            Figure size tuple of (width, height), in inches.

        :param padding:
            Amount of padding to place in between tracks, as a fraction of
            `yheight`
        """
        self.config = config
        self.figsize = figsize
        self.yheight = yheight
        self.padding = padding

        for features, kwargs in self.config:
            if "ybase" in kwargs:
                raise ValueError(
                    'Please do not specify "ybase"; this '
                    "is handled automatically by the %s class" % self.__class__.__name__
                )
            if "yheight" in kwargs:
                raise ValueError(
                    'Please do not specify "yheight", '
                    "this should be a separate arg to the %s "
                    "constructor" % self.__class__.__name__
                )

    def plot(self, ax=None):
        """
        If `ax` is None, create a new figure.  Otherwise, plot on `ax`.
        Iterates through the configuration, plotting each BedTool-like object
        as a separate track.
        """
        if ax is None:
            fig = plt.figure(figsize=self.figsize)
            ax = fig.add_subplot(111)
        yticks = []
        yticklabels = []
        ybase = 0
        i = 0
        padding = self.yheight * self.padding

        # Reverse config because incremental Track plotting works from bottom
        # up; this plots user-provided tracks in order from top down
        for features, kwargs in self.config[::-1]:
            t = Track(features, yheight=self.yheight, ybase=ybase, **kwargs)
            ybase = t.ymax + padding
            ax.add_collection(t)
            if "label" in kwargs:
                yticklabels.append(kwargs["label"])
            else:
                yticklabels.append(str(i))
                i += 1
            yticks.append(t.midpoint)

        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)

        ax.axis("tight")
        return ax


class BedToolsDemo(TrackCollection):
    def __init__(
        self,
        config,
        method,
        data_path=None,
        result_kwargs=None,
        method_kwargs=None,
        title_kwargs=None,
        new_style=True,
        subplots_adjust=None,
        *args,
        **kwargs
    ):
        """
        Class to handle BEDTools demos in a way that maintains flexibility.

        If the `config` list contains only one item, assume the method is one
        of the "-i" tools that only operate on one file.

        If the `config` list contains two items, then use the first as "-a" and
        the second as "-b".

        :param config:
            Either a list of (filename, options) tuples -- see docstring for
            TrackCollection for more info.

        :param method:
            Method of `BedTool` object to use, e.g., 'intersect'

        :param data_path:
            If not None, this path will be prepended to the files listed in
            `config`

        :param result_kwargs:
            Configuration for the results track.  This isn't added to the
            config list because the results haven't been created yet...

        :param method_kwargs:
            Keyword argument that are passed to the method, e.g., `u=True`

        :param title_kwargs:
            Keyword args for plot title (the text itself will come from the
            command that was run; this is for things like font size)

        :param new_style:
            Edit commands so that they use the "new style" BEDTools calls
            ("bedtools intersect" rather than "intersectBed")

        :param subplots_adjust:
            Additional kwargs sent to the figure's subplots_adjust() method,
            e.g., `dict(top=0.7)`


        :param args:
            Addtional arguments sent to TrackCollection

        :param kwargs:
            Additional keyword arguments sent to TrackCollection
        """
        if method_kwargs is None:
            method_kwargs = {}
        if result_kwargs is None:
            result_kwargs = {}
        if title_kwargs is None:
            title_kwargs = {}
        self.title_kwargs = title_kwargs
        self.new_style = new_style
        self.subplots_adjust = subplots_adjust

        # convert lists to tuples, cause we're going to edit the paths
        config = [list(i) for i in config]
        if data_path:
            for conf in config:
                if not isinstance(conf[0], six.string_types):
                    raise ValueError(
                        "data_path was specified, so you need "
                        "filenames in the config"
                    )
                conf[0] = os.path.join(data_path, conf[0])

        bt1 = pybedtools.BedTool(config[0][0])
        method = getattr(bt1, method)
        if len(config) == 2:
            result = method(config[1][0], **method_kwargs)
        elif len(config) == 1:
            result = method(**method_kwargs)
        else:
            raise ValueError(
                "`config` must have length 1 (for '-i' tools) or "
                "length 2 (for '-a -b' tools)."
            )

        config.append((result, result_kwargs))
        self.result = result
        super(BedToolsDemo, self).__init__(config, *args, **kwargs)

    def plot(self, ax=None):
        ax = super(BedToolsDemo, self).plot(ax)
        cmds = self.result._cmds[:]
        if self.new_style:
            cmds[0] = (
                "bedtools %s"
                % pybedtools.settings._prog_names[os.path.basename(cmds[0])]
            )
        ax.set_title(" ".join([os.path.basename(i) for i in cmds]), **self.title_kwargs)
        if self.subplots_adjust:
            ax.figure.subplots_adjust(**self.subplots_adjust)
        return ax


class ConfiguredBedToolsDemo(BedToolsDemo):
    def __init__(self, yaml_config, method, method_kwargs, **kwargs):
        """
        Wrapper around BedToolsDemo class that reads in a YAML config file.
        Useful for using the same "style" configuration many times.

        Contents of `yaml_config` must be YAML versions of BedToolsDemo args
        and kwargs **except** `method` and `method_kwargs`.
        """
        import yaml

        conf = yaml.load(open(yaml_config))

        disallowed = ["method", "method_kwargs"]
        for dis in disallowed:
            if dis in conf:
                raise ValueError("'%s' cannot be provided in the YAML config" % dis)

        conf["method"] = method
        conf["method_kwargs"] = method_kwargs
        conf.update(kwargs)
        super(ConfiguredBedToolsDemo, self).__init__(**conf)


if __name__ == "__main__":
    """
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
    """
    conf_file = pybedtools.example_filename("democonfig.yaml")
    data_path = pybedtools.example_filename("")  # dir name
    ax1 = ConfiguredBedToolsDemo(
        conf_file, method="intersect", method_kwargs={}, data_path=data_path
    ).plot()
    ax2 = ConfiguredBedToolsDemo(
        conf_file, method="intersect", method_kwargs=dict(u=True), data_path=data_path
    ).plot()
    plt.show()
