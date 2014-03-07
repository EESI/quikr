function xstar = quikrCustomTrained(trainingmatrix,inputfasta,k,lambda)
%xstar=quikrCustomTrained(traininmatrix,inputfasta,k,lambda) is the
%implementation of qukir that lets you use a custom training matrix
%"trainingmatrix" on the input data file "inputfasta" for the k-mer size
%"k" (note this k-mer size must be such that the number of rows of
%"trainingmatrix"=4^k), with the regularization value "lambda". The vector
%of predicted concentrations "xstar" is returned on the same basis as
%"trainingmatrix".
if nargin~=4
    error('There must be exactly 4 input arguments: the training matrix, the /path/to/input/fastafile, the k-mer size, and lambda');
end

[rows, columns]=size(trainingmatrix); %get the size of the training matrix
if rows~=4^k
    error('Wrong k-mer size for input training matrix');
end

[status, counts]=unix([sprintf('kmer_total_count -k %d', k) ' -i ' inputfasta]); %count the k-mers in the fasta file, in the forward direction, return the counts without labels.
if status ~= 0
  error('kmer_total_count failed: ensure kmer_total_count is in your path.');
end

counts=textscan(counts,'%f'); %read them in as floats.
counts=counts{:};
counts=counts/sum(counts); %normalize the counts into a probability vector
yaux=[0;lambda*counts]; %form the sample vector


Aaux=[ones(1,columns);lambda*trainingmatrix]; %form the k-mer sensing matrix
warning off
xstar=lsqnonneg(Aaux,yaux); %perform the non-negative lease squares
warning on
xstar=xstar/sum(xstar); %return the results as a probability vector 
