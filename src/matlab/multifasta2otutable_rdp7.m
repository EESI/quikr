% Quikr multifasta->otu_table_(for_qiime_use) wrapper code written by Gail Rosen -- 2/1/2013
%This is an example of how to run Quikr on the default RDP_7 training set

%make sure Matlab/Octave is in your path
%cd /path/to/Quikr

%User-defined variables
input_directory='../../separated_samples'; %path to input directory of samples
output_directory='quikr_results'; %path to where want output files to go
otu_table_name='rdp_otu_table.txt'; %name of output otu_table filename
%Do not have to define trainingdatabase file here
[headers,~]=fastaread('../../data/trainset7_112011.fa'); %read in the training database

mkdir([output_directory])
thedirs=dir([input_directory]);
thetime=zeros(numel(thedirs)-1,1);
names={};


if(exist('OCTAVE_VERSION')) %check to see if running Octave or Matlab

%This is Octave Version


species=struct();
keys={};

tic()

i=0;
%for numdirs=3:5
for numdirs=3:numel(thedirs)
i=i+1;
disp([num2str(i) ' out of ' num2str(numel(thedirs)-2)])
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
  names{j}=namesandproportions{1,j}(1:strfind(namesandproportions{1,j},' ')-1);
  names{j}=strrep(names{j},'|','_');
  amount(j)=namesandproportions{2,j};
if strcmp(keys,names{j})
	temp=species.(names{j});
      	temp(i)=round(amount(j).*numreads);
      	species.(names{j})=temp;
  else
      temp=zeros(numel(thedirs)-3+1,1);
      temp(i)=round(amount(j).*numreads);
      if temp(i)==0
	% do not make a key -- has insignificant counts
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

tic();
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

for k=1:numel(keys)
truncname=strrep(keys{k},'_','|');
  fprintf(fid,'%s',truncname);

  temp(:,k)=species.(keys{k});
        for i=1:numits
		fprintf(fid,'\t%d',temp(i,k));
	end
fprintf(fid,'\n');
end
fclose(fid);

disp('Time to output OTU Table:')
toc()

else

%This is the Matlab version


species=containers.Map;

tic()

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
      if temp(i)==0
	% insignificant counts -- do not put
      else
        species(names{j})=temp;
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

'Total time to compute Quikr:'
toc()
'Quikr Average time per file:'
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

'Time to output OTU Table:'
toc

end
