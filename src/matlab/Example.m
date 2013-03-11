% Here is an example of re-training the quikr method, performing the reconstruction, and outputing the results to a CSV file appropriate for Excel.
mat=quikrTrain('TaxCollectortrainset7_112011.fa',6); %retrain using the TaxCollector version of trainset7_112011.fa
xstar=quikrCustomTrained(mat,'testfastafile',6,10000); %perform the reconstruction
quikrWriteOutput('/path/to/output/director','filenameforoutput',full(xstar),'/path/to/TaxCollectortrainset7_112011.fa',2) %write the output in CSV format, the last argument (2) corresponds to the phylum level. See the file quikrWriteOutput.m for more details.




%This is an example of how to run Quikr

%If you want to use the default training-set (that is, trainset7_112011.fa
%from RPD), then execute the following

%make sure Matlab/Octave is in the Quikr director
cd /path/to/Quikr

fastafilename='/path/to/your/fastafile.fasta'; %full path name to your data file
xstar=quikr(fastafilename); %this will give the predicted reconstruction frequencies


%xstar will be on the same basis as trainset7_112011.fa, so to get the
%sequences that are predicted to be present in your sample:
[headers,~]=fastaread('trainset7_112011.fa'); %read in the training database
%note fastaread is not by default included in Octave, the fastaread.m file
%is included in the Quikr download directory however, and is directly
%compatible with Matlab.
nonzeroentries=find(xstar); %get the indicies of the sequences quikr predicts are in your sample
proportionscell=num2cell(xstar(nonzeroentries)); %convert the concentrations into a cell array
namescell=headers(nonzeroentries); %Get the names of the sequences
namesandproportions={namescell{:}; proportionscell{:}}; %This cell array contains the (unsorted) names of the reconstructed sequences and their concentrations (in the first and second columns respectively)
%so to find which sequence is the most abundant in your mixture:
[val,ind]=max(xstar(nonzeroentries)); %get the maximum value and it's position
namesandproportions{1:2,ind} %note that this does not imply this specific strain or species is in your sample, just that phylum/class/order/family/genus this species belongs to is in your sample.


% If you would like to use a custom training database, follow the following
% steps
fastafilename='/path/to/your/fastafile.fasta'; %full path name to your data file
trainingdatabasefilename='/path/to/your/trainingdatabase.fasta'; %full path to the FASTA file you wish to use as a training database
k=6; %pick a k-mer size
trainingmatrix=quikrTrain(trainingdatabasefilename,k); %this will return the training database
%then to do the reconstuction
lambda=10000; %pick a lambda (larger lambda -> theoretically predicted concentrations are closer to actual concentrations), this depends on k-mer size picked, also size and condition of the TrainingMatrix
xstar=quikrCustomTrained(trainingmatrix,fastafilename,k,lambda); %get the predicted reconstruction frequencies

%again Xstar is on the same basis as the TrainingMatrix, so to get the
%sequences that are predicted to be present in your sample:
[headers,~]=fastaread(trainingdatabasefilename); %read in the training database
nonzeroentries=find(xstar); %get the indicies of the sequences quikr predicts are in your sample
proportionscell=num2cell(xstar(nonzeroentries)); %convert the concentrations into a cell array
namescell=headers(nonzeroentries); %Get the names of the sequences
namesandproportions={namescell{:}; proportionscell{:}}; %This cell array contains the (unsorted) names of the reconstructed sequences and their concentrations (in the first and second columns respectively)
