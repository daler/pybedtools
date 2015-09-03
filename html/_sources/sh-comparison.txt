.. include:: includeme.rst

.. _shell_comparison:

Shell script comparison
=======================
The following two scripts illustrate the same analysis.  The first script uses
:mod:`pybedtools`, and the second uses bash scripting.  The filenames in these
scripts are written so they can be run without modification from the
`pybedtools/scripts` source directory.

Both scripts print the genes that are <5000 bp from intergenic SNPs.  These
scripts show how the same analysis can be performed with :mod:`pybedtools` in
a much clearer and reusable fashion without losing any speed. Furthermore, note
that the bash script requires knowledge in three languages (Perl, bash, and
awk) to accomplish the same thing as the Python script.

The `bash` script contains comparative notes as well as timing comparisons
between the two scripts.


:mod:`pybedtools` version (`py_ms_example.py`)
----------------------------------------------
.. literalinclude:: ../../pybedtools/scripts/py_ms_example.py


`bash` version (`sh_ms_example.sh`)
-----------------------------------
.. literalinclude:: ../../pybedtools/scripts/sh_ms_example.sh
    :language: bash

