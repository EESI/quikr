# Quikr's Matlab Implementation #
The Quikr implementation works in Matlab and also works in Octave, but the
Octave version will run much slower

## Quikr Example ##
This is an example of how to run Quikr. Before you try the example please make
make sure that you are in the quikr's matlab directory (src/matlab/):

    cd quikr/src/matlab

### Using Quikr with the default databse ###
This is the full path name to your data file:

    fastafilename='/path/to/quikr-code/testfastafile.fa';

This will give the predicted reconstruction frequencies using the default
training database trainset7\_112011.fa from RDP version 2.4 Xstar will be on the
same basis as trainset7\_112011.fa, so to get the sequences that are predicted
to be present in your sample:

    xstar=quikr(fastafilename);

Read in the training database.

_Note fastaread is not by default included in Octave. The fastaread.m file is
included in the Quikr download directory and is directly compatible with
Matlab._

    [headers,~]=fastaread('trainset7\_112011.fa');


Get the indicies of the sequences quikr predicts are in your sample

    nonzeroentries=find(xstar);

Convert the concentrations into a cell array

    proportionscell=num2cell(xstar(nonzeroentries));

Get the names of the sequences

    namescell=headers(nonzeroentries);

This cell array contains the (unsorted) names of the reconstructed sequences and
their concentrations (in the first and second columns respectively) so to find
which sequence is the most abundant in your mixture:

    namesandproportions={namescell{:}; proportionscell{:}};

Get the maximum value and it's position

    [val,ind]=max(xstar(nonzeroentries));

Note that this does not imply this specific strain or species is in your sample,
just that phylum/class/order/family/genus this species belongs to is in your
sample.

    namesandproportions{1:2,ind}

### Using Quikr With A Custom Trained Database ###
If you would like to use a custom training database, follow the following steps:

Full path name to your data file: 

    fastafilename='/path/to/your/fastafile.fasta';

Full path to the FASTA file you wish to use as a training database

    trainingdatabasefilename='/path/to/your/trainingdatabase.fasta';

Pick a k-mer size (typically 6 )

    k=6;

This will return the training database then to do the reconstruction

    trainingmatrix=quikrTrain(trainingdatabasefilename,k);

Pick a lambda (larger lambda -> theoretically predicted concentrations are
closer to actual concentrations), this depends on k-mer size picked, also size
and condition of the TrainingMatrix

    lambda=10000;

Get the predicted reconstruction frequencies
Again, xstar is on the same basis as the TrainingMatrix, so to get the sequences
that are predicted to be present in your sample:

    xstar=quikrCustomTrained(trainingmatrix,fastafilename,k,lambda);

Read in the training database

    [headers,~]=fastaread(trainingdatabasefilename);

Get the indices of the sequences quikr predicts are in your sample

    nonzeroentries=find(xstar);

Convert the concentrations into a cell array

    proportionscell=num2cell(xstar(nonzeroentries));

Get the names of the sequences

    namescell=headers(nonzeroentries);

This cell array contains the (unsorted) names of the reconstructed sequences and
their concentrations (in the first and second columns respectively)

    namesandproportions={namescell{:}; proportionscell{:}};
    
### Using Multifasta2otu.m ###

Usage tips:
* Please name fasta files of sample reads with <sample id>.fa<*> and place them
  into one directory without any other f ile in that directory (for example, no
  hidden files that the operating system may generate, are allowed in that
  direct ory)
* Note: When making your QIIME Metadata file, the sample id's must match the
  fasta file prefix names
* Fasta files of reads must have a suffix that starts with .fa (e.g.: .fasta and
  .fa are valid while .fna is NOT)
* Modify the top of the Matlab/Octave scripts for <input_directory>,
  <output_directory>, <output_filename>, and <train ing_database_filename>

To use with QIIME, one must run the QIIME conversion tool on our OTU table
output:

    convert_biom.py -i <quikr_otu_table.txt> -o <quikr_otu>.biom
    --biom_table_type="otu table"


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
