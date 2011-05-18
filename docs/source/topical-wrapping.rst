Wrapping new tools
==================
This section serves as a reference for wrapping new tools as they are added to
BEDTools.

Method name
-----------
Generally, I've tried to keep method names as similar as possible to
BEDTools programs while still being PEP8-compliant.  The trailing 'Bed' is
usually removed from the program name.

Args and kwargs
---------------
If a BEDTools program accepts a `-b` argument, then the signature and first
line should look like this::

    def new_method(self, b=None, **kwargs):
        kwargs['b'] = b
        ...

This is to allow another BedTool (or file, or stream) be passed as the
first non-keyword argument; otherwise, the user would always have to do::

    a.intersect(b=b)
 
instead of::

    a.intersect(b)

It's a minor thing, but it's convenient.  Furthermore, the line::

    kwargs['b'] = b

just makes sure that `b` is in the kwargs dict for when it's passed to the
:meth:`BedTool.handle_kwargs` method.

If the program you're wrapping doesn't accept another BED-like file as
`-b`, then the method should only accept `self` and `**kwargs`::

    def new_method(self, **kwargs):
        ...


.. TODO:: left off here....

Setting the default
-------------------
Generally, the next thing to do is to set the default, or implicit kwarg.
This is generally the `-i` or `-a` kwarg.


