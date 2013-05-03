%This is an example of how to run Quikr on the default RDP_7 training set, trainset7_112011.fa

%make sure Matlab/Octave is in your path
%cd /path/to/Quikr

%User-defined variables
input_directory='../../separated_samples'; %path to input directory of samples
output_directory='quikr_results'; %path to where want output files to go
otu_table_name='rdp_otu_table.txt'; %name of output otu_table filename
%Do not have to define trainingdatabasename here

mkdir([output_directory])
thedirs=dir([input_directory]);
thetime=zeros(numel(thedirs)-1,1);
names={};

species=containers.Map;

tic()

[headers,~]=fastaread('trainset7_112011.fa'); %read in the training database

i=0;
%for numdirs=3:6
for numdirs=3:numel(thedirs)
i=i+1;
[num2str(i) ' out of ' num2str(numel(thedirs)-2)]
fastafilename=[input_directory '/' thedirs(numdirs).name];
[loadfasta,~]=fastaread(fastafilename);
numreads=numel(loadfasta);
xstar=quikr(fastafilename); %this will give the predicted reconstruction frequencies

nonzeroentries=find(xstar); %get the indicies of the sequences quikr predicts are in your sample
proportionscell=num2cell(xstar(nonzeroentries)); %convert the concentrations into a cell array
namescell=headers(nonzeroentries); %Get the names of the sequences
namesandproportions={namescell{:}; proportionscell{:}}; %This cell array contains the (unsorted) names of the reconstructed sequences and their concentrations (in the first and second columns respectively)

[a cols]=size(namesandproportions);
amount=zeros(cols,1);
for j=1:cols
  names{j}=namesandproportions{1,j};
  amount(j)=namesandproportions{2,j};
  if isKey(species,names{j})
      temp=species(names{j});
      temp(i)=round(amount(j).*numreads);
      species(names{j})=temp;
  else
      temp=zeros(numel(thedirs)-3+1,1);
      temp(i)=round(amount(j).*numreads);
      species(names{j})=temp;
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

'Total time to compute Quikr'
toc()
'Quickr Average time per file'
mean(diff(thetime(1:i+1)))

tic
numits=i;


fid=fopen([output_directory '/' otu_table_name],'w');
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

thekeys=species.keys;
for k=1:species.Count
if k==1
  space=thekeys{k}(20);
end
delimit =strfind(thekeys{k},space);
truncname=thekeys{k}(1:delimit-1);
  fprintf(fid,'%s',truncname);

  temp(:,k)=species(thekeys{k});
        for i=1:numits
		fprintf(fid,'\t%d',temp(i,k));
	end
fprintf(fid,'\n');
end
fclose(fid);

'Time to output OTU Table'
toc
