function xstar = quikr(inputfasta)
%xstar=quikr('path/to/input/fasta/file') returns the estimated frequencies
%of bacteria present when given an input FASTA file of amplicon (454)
%reads. A k-mer based, L1-regularized, sparsity promoting algorithm is
%utilized. This is the default implementation of Quikr, here lambda=10,000, k-mer
%size is 6, and we use the trainset7_112011 database from RDP to estimate
%the present species. In practice, reconstruction is accurate only down to
%the genus level (not species or strain).
%The program ./count-linux or ./count-osx needs to be placed in the same 
%directory as this script. In Linux, you might need to: chmod a+rx count.
%Also the inputfasta needs to be a standard FASTA file. The input FASTA file must have reads all in the same
%orientation, so as a preprocessing step, make sure everything is in the
%forward (+) orientation.
if nargin>1
  error('too many input arguments');
end

[status, counts]=unix(['count-kmers -r 6 -1 -u ' inputfasta]); %count the 6-mers in the fasta file, in the forward direction, return the counts without labels
if status ~= 0
  error('count-kmers failed: ensure count-kmers is in your path.');
end

counts=textscan(counts,'%f'); %convert into floats
counts=counts{:}; %make into a vector
counts=counts/sum(counts); %form a probability vector from the counts
lambda=10000; %this is the default lambda value (see the paper), and is used in the formation of trainset7_112011N6Aaux.mat
yaux=[0;lambda*counts]; % form the sample vector
load('../../data/trainset7_112011N6Aaux.mat') %load the 6-mer sensing matrix...I need to rename this so the variable is just Aaux.
if not(exist('Aaux2')) %make sure Aaux2 loaded
    error('The default training database trainset7_112011N6Aaux.mat failed to load. Be sure it exists in the Quikr folder or current directory');
end
%based on RDP's trainset 7 from 11/2011
isOctave = exist('OCTAVE_VERSION') ~= 0; %check to see if running Octave or Matlab
if isOctave %if octave, need to manually read in sparse-matrix format
    n=length(Aaux2.jc)-1;
    ii=repelems(1:n,[1:n; diff(Aaux2.jc)]);
    Aaux=sparse(Aaux2.ir+1,ii,Aaux2.data);
else
    Aaux=Aaux2; %fix the names
    clear Aaux2;
end %Matlab automatically reads it in correctly, just need to rename it
warning off
xstar=lsqnonneg(Aaux,yaux); %perform the nonnegative least squares
warning on
xstar=xstar/sum(xstar); %transform the output into a probability vector. Note this vector is on the same basis as Aaux, so the entries correspond to sequences in trainset7_112011.fa
