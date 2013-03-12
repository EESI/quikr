# Quikr Command Line Utilities #  

Quikr has three command-line utilities that mirror the behavior of the python
module and the matlab implementation. The advantage of this is ease of scripting
and job management. These utilities are written in python and wrap the quikr
module.

## quikr\_train ##

The quikr\_train is a tool to train a database for use with the quikr tool.
Before running the quikr utility, you need to generate the trained matrix or
download a pretrained matrix from our database\_download.html.

### Usage ###
quikr\_train returns a custom trained matrix that can be used with the quikr
function. You must supply a kmer.

quikr\_train's optional arguments:
  -i, --input, the database of sequences (fasta format)
  -o, --output, the trained matrix (text file)
  -k, --kmer, the kmer size (integer)
  -z, --compress  compress the output matrix with gzip (flag)

## quikr ##
Quikr returns the estimated frequencies of batcteria present when given a
input FASTA file. A default trained matrix will be used if none is supplied
You must supply a kmer and default lambda if using a custom trained matrix.

quikr's optional arguments:
  -f, --fasta, the fasta file sample
  -o, --output OUTPUT, the output path (csv output)
  -t, --trained-matrix, the trained matrix
  -l, --lamb, the lambda size. (the default lambda value is 10,000)
  -k, --kmer, this specifies which kmer to use (default is 6)


### Troubleshooting ###

If you are having trouble, and these solutions don't work. Please contact the
developers with questions and issues.

#### Broken Pipe Errors #### 
Make sure that you have the count-kmers and probablilties-by-read in your
$PATH, and that they are executable. 

If you have not installed quikr system-wide, you'll need to add the folder
location of these binaries in the terminal before running the command:
 
    PATH = $PATH:/path/to/quikr/src/nbc/

Make sure that the binaries are executable by running:

    chmod +x probabilities-by-read
    chmod +x count-kmers
   
#### Python Cannot Find XYZ ####

Ensure that you have Python 2.7, Scipy, Numpy, and BIOpython installed 
and that python is setup correctly. You should be able to do this from a python
prompt without any errors:
    >>> import numpy
    >>> import scipy
    >>> from Bio import SeqIO

