# Quikr Command Line Utilities #
Quikr has three command-line utilities that mirror the behavior of the python
module and the matlab implementation. The advantage of this is ease of scripting
and job management, as well as faster processing and lower memory usage. These 
utilities are written in C and utilize OpenMP for multithreading.

For more in-depth information about these tools please refer to their
respective manual pages.

## Quikr\_train ##
The quikr\_train is a tool to train a database for use with the quikr tool.
Before running the quikr utility, you need to generate the sensing matrix.

### Usage ###
quikr\_train returns a custom sensing matrix that can be used with the quikr
function.

    quikr_train's arguments:
      -i, --input, the database of sequences (fasta format)
      -o, --output, the sensing matrix (text file)
      -k, --kmer, specifiy wha size of kmer to use. (default value is 6)
      -v, --verbose, verbose mode.
      -V, --version, print version.

### Example ###
Here is an example on how to train a database. This uses the -z flag to compress
the output matrix since it can be very large. Because of the sparse nature of
the database, the matrix easily achieves a high compression ratio, even with
gzip. It takes the gg94\_database.fasta as an input and outputs the sensing 
matrix as gg94\_sensing\_databse.npy.gz

    quikr_train -i gg94_database.fasta -o gg94_sensing_database.matrix.gz -k 6

## Quikr ##
Quikr returns the estimated frequencies of batcteria present when given a
input FASTA file. You need to train a matrix or download a new matrix 

### Usage ###
quikr returns the solution vector as a csv file.

    quikr's arguments:
    -i, --input the sample's fasta file of NGS READS (fasta format)
    -f, --sensing-fasta location of the fasta file database used to create the sensing matrix (fasta format)
    -s, --sensing-matrix location of the sensing matrix. (trained from quikr_train)
    -k, --kmer specify what size of kmer to use. (default value is 6)
    -l, --lambda lambda value to use. (default value is 10000)
    -o, --output OTU_FRACTION_PRESENT a vector representing the percentage of database sequence's presence in sample. (csv output)
    -v, --verbose verbose mode.
    -V, --version print version.
    -d, --debug debug mode, read manpage for more details

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

    multifasta_to_otu's arguments:
    -i, --input-directory the directory containing the samples' fasta files of 
    reads (note each file should correspond to a separate sample)
    -f, --sensing-fasta location of the fasta file database used to create the sensing matrix (fasta format)
    -s, --sensing-matrix location of the sensing matrix. (sensing from quikr_train)
    -k, --kmer specify what size of kmer to use. (default value is 6)
    -l, --lambda lambda value to use. (default value is 10000)
    -j, --jobs specifies how many jobs to run at once. (default value is the number of CPUs)
    -o, --output the OTU table, with NUM_READS_PRESENT for each sample which 
    is compatible with QIIME's convert_biom.py (or a sequence table if not OTU's)
    -v, --verbose verbose mode.
    -V, --version  print version.

### Post-processing of Multifasta\_to\_otu  ###

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


#### Broken Pipe Errors ####
Make sure that you have the count-kmers and probablilties-by-read in your
$PATH, and that they are executable. 

If you have not installed quikr system-wide, you'll need to add the folder
location of these binaries in the terminal before running the command:
 
    mv /path/to/quikr/src/nbc/count /path/to/quikr/src/nbc/count-kmers
    PATH = $PATH:/path/to/quikr/src/nbc/
