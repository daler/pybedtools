import os
import sys
import sqlite3
import pybedtools
import time
import collections


def now():
    return time.time()


def get_name(fname):
    return os.path.splitext(os.path.basename(fname))[0]


class IntersectionMatrix(object):
    """
    Class to handle many pairwise comparisons of interval files
    """

    def __init__(self, beds, genome, iterations, dbfn=None, force=False):
        """
        Class to handle and keep track of many pairwise comparisons of interval
        files.

        A lightweight database approach is used to minimize computational time.

        The database stores filenames and calculation timestamps;
        re-calculating a matrix using the same interval files will only
        re-calculate values for those files whose modification times are newer
        than the timestamp in the database.

        `beds` is a list of bed files.

        `genome` is the string assembly name, e.g., "hg19" or "dm3".

        `dbfn` is the filename of the database you'd like to use to track
        what's been completed.

        Example usage:

        First, get a list of bed files to use:
        #>>> beds = [
        #... pybedtools.example_filename(i) for i in  [
        #... 'Cp190_Kc_Bushey_2009.bed',
        #... 'CTCF_Kc_Bushey_2009.bed',
        #... 'SuHw_Kc_Bushey_2009.bed',
        #... 'BEAF_Kc_Bushey_2009.bed'
        #... ]]

        Set some parameters.  "dm3" is the genome to use; info will be stored
        in "ex.db".  `force=True` means to overwrite what's in the database
        #>>> # In practice, you'll want many more iterations...
        #>>> im = IntersectionMatrix(beds, 'dm3',
        #...            dbfn='ex.db', iterations=3, force=True)
        #>>> # Use 4 CPUs for randomization
        #>>> matrix = im.create_matrix(verbose=True, processes=4)
        """
        self.beds = beds
        self.genome = genome
        self.dbfn = dbfn
        self.iterations = iterations

        if self.dbfn:
            self._init_db(force)
            self.conn = sqlite3.connect(dbfn)
            self.conn.row_factory = sqlite3.Row
            self.c = self.conn.cursor()

    def _init_db(self, force=False):
        """
        Prepare the database if it doesn't already exist
        """
        if self.dbfn is None:
            return
        if os.path.exists(self.dbfn) and not force:
            return
        conn = sqlite3.connect(self.dbfn)
        c = conn.cursor()
        if force:
            c.execute("DROP TABLE IF EXISTS intersections;")
        c.executescript(
            """
        CREATE TABLE intersections (
            filea TEXT,
            fileb TEXT,
            timestamp FLOAT,
            actual FLOAT,
            median FLOAT,
            iterations INT,
            self INT,
            other INT,
            fractionabove FLOAT,
            fractionbelow FLOAT,
            percentile FLOAT,
            PRIMARY KEY (filea, fileb, iterations));
        """
        )
        conn.commit()

    def get_row(self, fa, fb, iterations):
        """
        Return the sqlite3.Row from the database corresponding to files `fa`
        and `fb`; returns None if not found.
        """
        if self.dbfn is None:
            return

        results = list(
            self.c.execute(
                """
                SELECT * FROM intersections
                WHERE
                filea=:fa AND fileb=:fb AND iterations=:iterations
                """,
                locals(),
            )
        )
        if len(results) == 0:
            return
        assert len(results) == 1
        return results[0]

    def done(self, fa, fb, iterations):
        """
        Retrieves row from db and only returns True if there's something in
        there and the timestamp is newer than the input files.
        """
        row = self.get_row(fa, fb, iterations)
        if row:
            tfa = os.path.getmtime(fa)
            tfb = os.path.getmtime(fb)
            if (row["timestamp"] > tfa) and (row["timestamp"] > tfb):
                return True
        return False

    def run_and_insert(self, fa, fb, **kwargs):
        a = pybedtools.BedTool(fa).set_chromsizes(self.genome)
        kwargs["iterations"] = self.iterations
        results = a.randomstats(fb, **kwargs)
        self.add_row(results)

    def add_row(self, results):
        """
        Inserts data into db.  `results` is a dictionary as returned by
        BedTool.randomstats with keys like::

            'iterations'
            'actual'
            'file_a'
            'file_b'
            self.fn
            other.fn
            'self'
            'other'
            'frac randomized above actual'
            'frac randomized below actual'
            'median randomized'
            'normalized'
            'percentile'
            'lower_%sth' % lower_thresh
            'upper_%sth' % upper_thresh
        """
        # translate results keys into db-friendly versions
        translations = [
            ("file_a", "filea"),
            ("file_b", "fileb"),
            ("median randomized", "median"),
            ("frac randomized above actual", "fractionabove"),
            ("frac randomized below actual", "fractionbelow"),
        ]
        for orig, new in translations:
            results[new] = results[orig]

        results["timestamp"] = now()

        sql = """
        INSERT OR REPLACE INTO intersections (

            filea,
            fileb,
            timestamp,
            actual,
            median,
            iterations,
            self,
            other,
            fractionabove,
            fractionbelow,
            percentile)

            VALUES (

            :filea,
            :fileb,
            :timestamp,
            :actual,
            :median,
            :iterations,
            :self,
            :other,
            :fractionabove,
            :fractionbelow,
            :percentile)

        """
        self.c.execute(sql, results)
        self.conn.commit()

    def create_matrix(self, verbose=False, **kwargs):
        """
        Matrix (implemented as a dictionary), where the final values are
        sqlite3.ROW objects from the database::

            {
                filea: {
                            filea: ROW,
                            fileb: ROW,
                            ...},
                fileb: {
                            filea: ROW,
                            fileb: ROW,
                            ...},

                        }
            }
        """
        nfiles = len(self.beds)
        total = nfiles ** 2
        i = 0
        matrix = collections.defaultdict(dict)
        for fa in self.beds:
            for fb in self.beds:
                i += 1

                if verbose:
                    sys.stderr.write("%(i)s of %(total)s: %(fa)s + %(fb)s\n" % locals())
                    sys.stderr.flush()

                if not self.done(fa, fb, self.iterations):
                    self.run_and_insert(fa, fb, **kwargs)

                matrix[get_name(fa)][get_name(fb)] = self.get_row(
                    fa, fb, self.iterations
                )

        return matrix

    def print_matrix(self, matrix, key):
        """
        Prints a pairwise matrix of values. `matrix` is a dict-of-dicts from
        create_matrix(), and `key` is a field name from the database -- one of:

        ['filea', 'fileb', 'timestamp', 'actual', 'median', 'iterations',
        'self', 'other', 'fractionabove', 'fractionbelow', 'percentile']
        """
