# How Can I Install The Quikr Utility? #
To use Quikr there are several prerequisites. 

## Requirements ##
+ Mac OS X 10.6.8 or GNU/Linux
+ 4Gb of RAM minimum. Absolutely neccessary.
+ gcc that supports OpenMP

### Python Requirements ###
*the python version is out of date*
+ Python 2.7+ 
+ Scipy 
+ Numpy
+ BioPython

### Matlab Requrements ###
+ [dna\_utils](http://github.com/EESI/dna-utils/) must be installed

### Mac Requirements ###
+ Mac OS X 10.6.8 (what we have tested)
+ GCC 4.7 or newer. (gcc 4.2 did not work, and is the default installation)
+ OpenMP libraries (libgomp, usually comes with gcc)

### Linux Requirements ###
+ GCC 4.7 or newer
+ OpenMP libraries (libgomp, usually comes with gcc)

We also have a Quikr implementation in Matlab so that you can easily integrate
Quikr into your custom programs and scripts.

## Installation ##
Our Quikr code is available on our sourceforge download page:

[sourceforge project page](http://sourceforge.net/projects/quikr/)

Alternatively you can clone the latest version from Github, or 
download older versions using github's 'release' view.

To install quikr, download our project and in the folder run:

    make
    sudo make install

This will install the quikr, quikr\_train and multifasta\_to\_otu utilities.
To install the python scripts and module systemwide, run

    make python
    sudo make install_python
