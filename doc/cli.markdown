# Quikr Command Line Utilities #  

Quikr has three command-line utilities that mirror the behavior of the python
module and the matlab implementation. The advantage of this is ease of scripting
and job management. These utilities are written in python and wrap the quikr
module.

## Quikr\_train ##

The quikr\_train is a tool to train a database for use with the quikr tool.
Before running the quikr utility, you need to generate the sensing matrix or
download a pretrained matrix from our database\_download.html.

### Usage ###
quikr\_train returns a custom trained matrix that can be used with the quikr
function. You must supply a kmer.

quikr\_train's arguments:
  -i, --input, the database of sequences (fasta format)
  -o, --output, the trained matrix (text file)
  -k, --kmer, the kmer size, the default is 6 (integer)
  -z, --compress  compress the output matrix with gzip (flag)

### Example ###
Here is an example on how to train a database. This uses the -z flag to compress
the output matrix since it can be very large. Because of the sparse nature of
the database, the matrix easily achieves a high compression ratio, even with
gzip. It takes the gg94\_database.fasta as an input and outputs the trained 
matrix as gg94\_trained\_databse.npy.gz

    quikr_train -i gg94_database.fasta -o gg94_trained_database.npy.gz -k 6 -z 

## Quikr ##
Quikr returns the estimated frequencies of batcteria present when given a
input FASTA file. A default trained matrix will be used if none is supplied
You must supply a kmer and default lambda if using a custom trained matrix.

### Usage ###
quikr returns the solution vector as a csv file.

quikr's arguments:
  -f, --fasta, the sample's fasta file of NGS READS
  -o, --output OTU\_FRACTION\_PRESENT, a vector representing the percentage of
  database sequence's presence in sample (csv output)
  -t, --trained-matrix, the trained matrix
  -l, --lamb, the lambda size. (the default lambda value is 10,000)
  -k, --kmer, this specifies the size of the kmer to use (default is 6)

## Multifasta\_to\_otu ##
The Multifasta\_to\_otu tool is a handy wrapper for quikr which lets the user
to input as many fasta files as they like, and then returns an OTU table of the
number of times a specimen was seen in all of the samples 

Warning: this program will use a large amount of memory, and CPU time. You can
reduce the number of cores used, and thus memory, by specifying the -j flag
with aspecified number of jobs. Otherwise python with run one job per cpu core.

# Pre-processing of Multifasta\_to\_otu  #

* Please name fasta files of sample reads with <sample id>.fa<*> and place them
  into one directory without any other file in that directory (for example, no
  hidden files that the operating system may generate, are allowed in that
  directory)
* Fasta files of reads must have a suffix that starts with .fa (e.g.: .fasta and
  .fa are valid while .fna is NOT)

### Usage ###
multifasta\_to\_otu's arguments:
  -i, --input-directory, the directory containing the samples' fasta files of
  reads (note each fasta file should correspond to a separate sample)
  -o, --otu-table, the OTU table, with OTU\_FRACTION\_PRESENT for each sample,
  which is compatible with QIIME's convert\_biom.py (or sequence table if not
  OTU's)
  -t, --trained-matrix, the trained matrix
  -f, --trained-fasta, the fasta file database of sequences
  -d, --output-directory, quikr output directory
  -l, --lamb, specify what size of lambda to use (the default value is 10,000)
  -k, --kmer, specify what size of kmer to use, (default value is 6)
  -j, --jobs, specifies how many jobs to run at once, (default=number of CPUs)

# Post-processing of Multifasta\_to\_otu  #

* Note: When making your QIIME Metadata file, the sample id's must match the
  sample fasta file prefix names

4-step QIIME procedure after using Quikr to obtain 3D PCoA graphs:
(Note: Our code works much better with WEIGHTED Unifrac as opposed to
Unweighted.)

Pre-requisites:
1. <quikr_otu_table.txt>
2. the tree of the database sequences that were used (e.g.  dp7\_mafft.fasttree,
   gg\_94\_otus\_4feb2011.tre, etc.)
3. your-defined <qiime_metadata_file.txt>

The QIIME procedue:
    convert_biom.py -i <quikr_otu_table.txt> -o <quikr_otu>.biom --biom_table_type="otu table"
    beta_diversity.py -i <quikr_otu>.biom -m weighted_unifrac -o beta_div -t <tree file> (example: rdp7_mafft.fasttree)>
    principal_coordinates.py -i beta_div/weighted_unifrac_<quikr_otu>.txt -o <quikr_otu_project_name>_weighted.txt
    make_3d_plots.py -i <quikr_otu_project_name>_weighted.txt -o <3d_pcoa_plotdirectory> -m <qiime_metadata_file>


# Python Quikr Troubleshooting #

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

