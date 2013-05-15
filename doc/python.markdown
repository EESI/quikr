# Python Documentation #
The python version comes with scripts that can be used like the regular quikr
program, and also a module called quikr so integration with python scripts 
is easier.

If you  are switching to use the python scripts instead of the regular, you
will need to regenerate your trained databases with the python version of 
quikr\_train.

## Function documentation ##
If you want to use our quikr module, run help on the module:

    >>> import quikr
    >>> help(quikr)
## Python Cannot Find XYZ ##

Ensure that you have Python 2.7, Scipy, Numpy, and BIOpython installed 
and that python is setup correctly. You should be able to do this from a python
prompt without any errors:
    >>> import numpy
    >>> import scipy
    >>> from Bio import SeqIO
