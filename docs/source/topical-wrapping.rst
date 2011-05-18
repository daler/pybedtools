Wrapping new tools
==================
This section serves as a reference for wrapping new tools as they are added to
BEDTools.

Let's assume we would like to wrap a new program, appropriately named
`newProgramBed`.  Its signature from the command line is `newProgramBed -a
<infile> -b <other file> [options]`, and it accepts `-a stdin` to indicate
data is being piped to it.

Method name
-----------
Generally, I've tried to keep method names as similar as possible to
BEDTools programs while still being PEP8-compliant.  The trailing 'Bed' is
usually removed from the program name. So here the name would probably be
`new_program`.  

Method signature (args and kwargs)
----------------------------------
If a BEDTools program accepts a `-b` argument (which is the case for this
example) then the signature and first line should look like this::

    def new_program(self, b=None, **kwargs):
        kwargs['b'] = b
        ...

That is, we specify a kwarg for `b` in the signature so that we have the
option calling this method either with a regular arg or a kwarg, and then
we ensure that `b` makes it into the `**kwargs` dictionary that will be
passed around later.

This is to allow another BedTool (or file, or stream) be passed as the
first non-keyword argument; otherwise, the user would always have to do::

    a.intersect(b=b)
 
instead of::

    a.intersect(b)

It's a minor thing, but it's convenient.

If the program you're wrapping doesn't accept another BED-like file as
`-b`, then the method should only accept `self` and `**kwargs`::

    def new_program(self, **kwargs):
        pass

Setting the default
-------------------
Generally, the next thing to do is to set the default, or implicit kwarg.
This is generally the `-i` or `-a` kwarg.

This is done by checking to see if the implicit kwarg was specified already
(which allows the user to override the implicit assumption) and if not, add
it to the kwargs dict.

For example::

    if 'a' not in kwargs:
        kwargs['a'] = self.fn

For BEDTools programs that optionally take an `-abam` argument, we need to
check to make sure that wasn't specified either::

    if ('abam' not in kwargs) and ('a' not in kwargs):
            kwargs['a'] = self.fn

So to continue the example, our method now looks like this::

    def new_program(self, b=None, **kwargs):

        if 'a' not in kwargs:
            kwargs['a'] = self.fn

Allowing input streams
----------------------
BEDTools programs that can accept `stdin` as their first input need to be
registered in the :meth:`BedTool.handle_kwargs` method, in the
`implicit_instream1` dictionary.

This dictionary specifies which keyword argument to look for an object that
could be used as the first input -- this could be a filename, an iterator
of features, or an open file.  Depending on what it finds there, it will
call the BEDTools program appropriately.  

For example, if :meth:`BedTool.handle_kwargs` finds an open file in
`kwargs['a']` for `intersectBed`, then it will tell `intersectBed` that
`-a` should be `stdin` and the open file will eventually be passed to the
`subprocess.Popen` call.  But if `kwargs['a']` is a filename, then it will
just tell `intersectBed` that `-a` should be the filename.

If there is a second argument that has the potential to be a stream (like
the `b` kwarg for :meth:`BedTool.intersect`, which can be a BedTool,
filename, IntervalIterator, or stream), then this kwarg should be added to
the `implicit_instream2` dictionary.

This second dictionary specifies which kwargs to check to see if they
contain an iterable.  If so, then it "renders" the iterable to a temp file
and will pass that filename to the BEDTools program. 

In the example case of :meth:`BedTool.intersect`, if `b` is a stream or an
iterable, the `handle_kwargs` method will convert that stream or iterable
to a tempfile (say, `tmp001`) and then will eventually send `intersectBed`
the option `-b tmp001`.

The :meth:`BedTool.handle_kwargs` method also decides whether to create a
new temp file (if `stream=False`) or not, and also converts kwargs like `a`
to `-a`.  Finally, it creates a list of commands ready for
:func:`call_bedtools` to run.

To illustrate, we add `a` to the `implicit_instream1` dict:

::

    implicit_instream1 = {'intersectBed': 'a',
                           'subtractBed': 'a',
                            'closestBed': 'a',
                             'windowBed': 'a',
                               'slopBed': 'i',
                              'mergeBed': 'i',
                               'sortBed': 'i',
                           'bed12ToBed6': 'i',
                            'shuffleBed': 'i',
                           'annotateBed': 'i',
                              'flankBed': 'i',
                          'fastaFromBed': 'bed',
                      'maskFastaFromBed': 'bed',
                           'coverageBed': 'a',
                           'newthingBed': 'a', # here's the new wrapped program
                           }

And back in the method body, we call :meth:`BedTool.handle_kwargs` so that
our method now looks like this::

    def new_program(self, b=None, **kwargs):

        if 'a' not in kwargs:
            kwargs['a'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='newthingBed', **kwargs)

Call BEDTools, and return result
--------------------------------
So now we have `cmds` (the list of commands ready to be called by
`subprocess.Popen`), `tmp` (the new tempfile ready to accept results, or
None if we will be streaming the output), and `stdin` (None if it's a
filename, or the open file that will be send to `subprocess.stdin`).  These
are passed to the :func:`call_bedtools` function, which does all the
`subprocess` business and returns a "stream", which could be any of the
things that a :class:`BedTool` can accept (filename, open file,
IntervalIterator).  Finally, we return a new :class:`BedTool` made from
this stream.

Now our method looks like this::

    def new_program(self, b=None, **kwargs):

        if 'a' not in kwargs:
            kwargs['a'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='newthingBed', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

Add decorators
--------------
Some decorators are used to add text to the method's docstring, like
`_implicit` and `_file_or_bedtool`, and `_returns_bedtool`.  More useful is
`_help`, which adds the full text of the BEDTools program to the end of the
docstring, so all information is available from the Python interpreter
(especially useful when using IPython).  The `_log_to_history` decorator
will register the calling of this method in the BedTool's history.

The final wrapped method, with all the decorators to add relevant text to
the docstring, then is simply::

    @_returns_bedtool()
    @_file_or_bedtool()
    @_help('newProgramBed')
    @_log_to_history
    @_implicit('-a')
    def new_program(self, b=None, **kwargs):

        if 'a' not in kwargs:
            kwargs['a'] = self.fn

        cmds, tmp, stdin = self.handle_kwargs(prog='newthingBed', **kwargs)
        stream = call_bedtools(cmds, tmp, stdin=stdin)
        return BedTool(stream)

Write tests!
------------
The only way to know for sure if your new wrapped method works is to write
good tests for it.  This can either be done in the docstring or in the test
suite.  See the source for how this is done; also check out the `test.sh`
script in the top level of the repository.

Send a pull request
-------------------
If you've made something that you think would be useful to others, please
send a github pull request so that your newly created and tested code can
be distributed to others.  Any contributions are much appreciated.
