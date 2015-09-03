Under the hood
==============

This section documents some details about what happens when a :class:`BedTool`
object is created and exactly what happens when a BEDTools command is called.
It's mostly useful for developers or for debugging.


There are three kinds of sources/sinks for BedTool objects:

* filename
* open file object
* iterator of Interval objects


Iterator "protocol"
-------------------
BedTool objects yield an Interval object on each `next()` call.  Where this
Interval comes from depends on how the BedTool was created and what format the
underlying data are in, as follows.

Filename-based
~~~~~~~~~~~~~~
If BED/GTF/GFF/VCF format, then use an `IntervalFile` object for Cython/C++
speed.

If SAM format, then use an `IntervalIterator`.  This is a Cython object that
reads individual lines and passes them to `create_interval_from_list`, a Cython
function.  `create_interval_from_list` does a lot of the work to figure out
what format the line is, and this is how we are able to support SAM Interval
objects.

If BAM format, then first do a Popen call to `samtools view`, and create an
`IntervalIterator` from subprocess.PIPE similar to SAM format.

Open file-based
~~~~~~~~~~~~~~~
All formats are passed to an `IntervalIterator`, which reads one line at
a time and yields an `Interval` object.

If it's a BAM file (specifically, a detected bgzip stream), then it's actually
first sent to the stdin of a `samtools` Popen call, and then the
subprocess.PIPE from that Popen's stdout is sent to an `IntervalIterator`.

Iterator or generator-based
~~~~~~~~~~~~~~~~~~~~~~~~~~~
If it's neither of the above, then the assumption is that it's already an
iterable of `Interval` objects.  This is the case if a `BedTool` is created
with something like::

    a = pybedtools.example_bedtool('a.bed')
    b = pybedtools.BedTool((i for i in a))


In this case, the `(i for i in a)` creates a generator of intervals from an
`IntervalFile` -- since `a` is a filename-based BedTool.  Since the first
argument to the BedTool constructor is neither a filename nor an open file, the
new BedTool `b`'s `.fn` attribute is directly set to this generator . . . so we
have a generator-based BedTool.

Calling BEDTools programs
-------------------------
Depending on the type of BedTool (filename, open file, or iterator), the method
of calling BEDTools programs differs.

In all cases, BEDTools commands are called via a `subprocess.Popen` call
(hereafter called "the Popen" for convenience).  Depending on the type of
BedTool objects being operated on, the Popen will be passed different objects
as stdin and/or stdout.

In general, using a filename as input is the most straightforward -- nothing is
passed to the Popen's stdin because the filenames are embedded in the BEDTools
command.

Using non-filename-based BedTools means that they are passed, one line at
a time, to the stdin of the Popen.  The commands for the BEDTools call
will specify "stdin" in these cases, as is standard for the BEDTools suite.

The default is for the output to be file-based.  In this case, an open tempfile
object is provided as the Popen's stdout.

If the returned BedTool is requested to be a "streaming" BedTool, then the
Popen's stdout will be subprocess.PIPE, and the new BedTool object will be
open-file based (which is what subprocess.PIPE acts like).

Specifically, here is the information flow of stdin/stdout for various
interconversions of BedTool types . . . .


:filename -> filename:
    The calling BedTool is filename-based and `stream=False`.

    * `stdin`: `None` (the filenames are provided in the BEDTools command)
    * `stdout`: open tempfile object
    * new BedTool: filename-based BedTool pointing to the tempfile's filename

:filename -> open file object:
    The calling BedTool is filename-based and `stream=True` is requested.

    * `stdin`: None (provided in the cmds)
    * `stdout`: open file object -- specifically, subprocess.PIPE
    * new BedTool: iterator-based BedTool.  Each `next()` call retrieves the
      next line in subprocess.PIPE

:open file object -> filename:
    The calling BedTool is from, e.g., subprocess.PIPE and there's
    a saveas() call to "render" to file.

    * `stdin`: each line in the open file object is written to subprocess.PIPE
    * `stdout`: open file object -- either a tempfile or new file created from
      supplied filename
    * new BedTool: filename-based BedTool

:open file object -> iterator:
    The calling BedTool is usually based on subprocess.PIPE, and the output
    will *also* come from subprocess.PIPE.

    * `stdin`: each line from the open file is written to subprocess.PIPE
    * `stdout`: open file object, subprocess.PIPE
    * new BedTool: filename based on subprocess.PIPE
