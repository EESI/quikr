% Quikr multifasta->otu_table_(for_qiime_use) wrapper code written by Gail Rosen -- 2/1/2013
%This is an example of how to run Multifasta Quikr with a custom 
%training database (in this case Greengenes OTU's within 94% identity)

%make sure Matlab/Octave is in your path
%cd /path/to/Quikr

%User-defined variables
input_directory='py-quikr/input/'; %path to input directory of samples
output_directory='quikr_results'; %path to where want output files to go
otu_table_name='gg1194_otu_octave.txt'; %name of output otu_table filename
trainingdatabasefilename='gg_94_otus_4feb2011.fasta'; %full path to the FASTA file you wish to use as a training database
k=6; %pick a k-mer size
addpath('quikr-code');

% make output directory
mkdir([output_directory]);

% get a list of the input directory

thedirs=dir([input_directory]);
thetime=zeros(numel(thedirs)-1,1);
names={};

tic();

trainingmatrix=quikrTrain(trainingdatabasefilename,k); %this will return the training database
disp('Training time:')
[headers,~]=fastaread(trainingdatabasefilename); %read in the training database
lambda=10000; 
training_time=toc()

species=struct();
keys={};

tic();


i=0;
%for numdirs=3:5
for numdirs=3:numel(thedirs)
i=i+1;
disp([num2str(i) ' out of ' num2str(numel(thedirs)-2)])
fastafilename=[input_directory '/' thedirs(numdirs).name];
[loadfasta,~]=fastaread(fastafilename);
numreads=numel(loadfasta);
xstar=quikrCustomTrained(trainingmatrix,fastafilename,k,lambda);

nonzeroentries=find(xstar); %get the indicies of the sequences quikr predicts are in your sample
proportionscell=num2cell(xstar(nonzeroentries)); %convert the concentrations into a cell array
namescell=headers(nonzeroentries); %Get the names of the sequences
namesandproportions={namescell{:}; proportionscell{:}}; %This cell array contains the (unsorted) names of the reconstructed sequences and their concentrations (in the first and second columns respectively)

[a cols]=size(namesandproportions);
amount=zeros(cols,1);
for j=1:cols
  names{j}=['s' namesandproportions{1,j}];
  amount(j)=namesandproportions{2,j};
  if strcmp(keys,names{j})
	temp=species.(names{j});
      	temp(i)=round(amount(j).*numreads);
      	species.(names{j})=temp;
  else
      temp=zeros(numel(thedirs)-3+1,1);
      temp(i)=round(amount(j).*numreads);
      if temp(i)==0
         %insignificant counts, do nothing
      else
         species.(names{j})=temp;
         keys{end+1}=names{j};
      end
  end
end

thefa=strfind(thedirs(numdirs).name,'.fa');

if ~isempty(thedirs(numdirs).name(1:thefa-1))
	sampleid{i}=thedirs(numdirs).name(1:thefa-1);
else
	sampleid{i}='empty_sampleid';
end

thetime(i+1)=toc();
thetime(i+1)

end

disp('Total time to compute Quikr:')
toc()
disp('Quikr Average time per file:')
mean(diff(thetime(1:i+1)))

tic()
numits=i;

fid=fopen([output_directory '/' otu_table_name],'w');

% insert our header
fprintf(fid,'# QIIME vGail OTU table\n');
fprintf(fid,'#OTU_ID\t');
for i=1:numits
  if i<numits
    fprintf(fid,'%s\t',sampleid{i});
  else
    fprintf(fid,'%s',sampleid{i});
  end
end

fprintf(fid,'\n');
%


% print our key and our values for each keys
for k=1:numel(keys)
  fprintf(fid,'%s',keys{k}(2:end))
  temp(:,k)=species.(keys{k});
 
  for i=1:numits
    fprintf(fid,'\t%d',temp(i,k));
  end
  
  fprintf(fid,'\n');
end

fclose(fid);

disp('Time to output OTU Table:')
toc()

