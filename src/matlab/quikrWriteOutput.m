function quikrWriteOutput(directory,filename,xstar,pathtotrainingdatabase,taxonomicrank)
%taxonomic rank: 1=kindom, 2=phylum, 3=class, 4=order, 5=family,
%6=genus, 7=species, 8=strain
%quikrWriteOutput(directory,filename,xstar,trainingdatabase,taxonomicrank)
%will output a Comma Seperated (CSV) file with the file name 'filename' in the directory
%'directory' that gives the summary of concentrations of the reconstruction
%at the taxonomic rank specified by 'taxonomicrank'. The database located 
%at'pathtotrainingdatabase' must be in the fasta format with the headers
%in the format output by TaxCollector.
%(see Giongo et al, 2010 TaxCollector: Modifying Current 16S rRNA
%Databases for the Rapid Classification at Six Taxonomic Levels, Diversity.
%doi:10.3390/d2071015).
%The output excel file is suitable for use in visualiztion software such as
%Krona (Ondov, Bergman, and Phillippy, 2011, Interactive Metagenomic
%Visualiztion in a Web Browser, BMC bioinformatics,
%doi:10.1186/1471-2105-12-385).
if nargin>5
    error('too many input arguments. Needs: directory, output filename, xstar, path to training database, and taxonomic rank');
end
if nargin<5
    error('too few input arguments. Needs: directory, output filename, xstar, path to training database, and taxonomic rank');
end
[headers,~]=fastaread(pathtotrainingdatabase);%read in the headers of the training database

for i=1:length(headers)
    if not(length(regexp(headers{2},'\[\d\]'))==9)
        error(sprintf('The header in entry %d in the training database is not in the Tax Collector format: eg ">[1]kingdom[2]phylum[3]class[4]order[5]family[6]genus[7]species[8]strain|accession"',i));
    end
end


seqnames=headers(find(xstar)); %get the names of the sequences that are given nonzero concentration in the reconstruction xstar
basis=find(xstar); %the indicies that have nonzero concentrations
        %seqnamesnew=cell(length(basis),1);
for i=1:length(basis)
    seqnamesnew{i}=strcat(sprintf('%f[0] ',xstar(basis(i))),seqnames{i}); %format is concentration, and then the name as it is in the training database
end
        
for i=1:length(seqnamesnew)
    splitup(i,:)=regexp(seqnamesnew{i},'\[\d\]','split'); %split everything up on the token [%d] as per the TaxCollector format
end
%xlswrite(fullfile(directory,[filename '.xls']),splitup);
%so there's a problem: trainset7_112011 is not in the tax collector format.
%So one thing I could do is use the accessions or something like that to
%get the corresponding headers from TaxCollectorRDP10_28 and create a new
%trainset7_112011 (or at least, one with new headers that then work
%according to TaxCollector).

%So it does look like that for each header in trainset7_112011, after the >
%and before the | that "accession" also appears at the end of some header in
%TaxCollectorRDP10_28 (after the | on the end). So I'd like to make a new
%trainset7_112011 in which the headers are simply replaced.


%Now I'm going to try to unique at the given taxanomic rank, and then tally
%up the concentrations of all the corresponding entries in splitup
ranknames=unique(splitup(:,taxonomicrank+2)); %get the unique names at the given taxonomic rank
concentrations=str2double(splitup(:,1)); %pull off the concentrations
       
for i=1:length(ranknames) %for each of the names
    selectlist=strcmp(splitup(:,taxonomicrank+2),ranknames(i)); %get the rows that have the given name in the proper taxonomic rank
    %newsplitup(i,:)={num2str(sum(concentrations(selectlist)),'%f'),ranknames(i)}; %sum up the concentrations in all these rows and put that next to the names just down to the given taxonomic rank
    newsplitup(i,:)={num2str(sum(concentrations(selectlist)),'%f'),splitup{find(selectlist,1),3:taxonomicrank+2}};
end
%xlswrite(fullfile(directory,[filename '.xls']),newsplitup); %write this as the final output

fid=fopen(fullfile(directory,[filename '.csv']),'w');
if fid==-1
    error(sprintf('The file %s is currently write-protected (e.g. open or in use in a different program)', fullfile(directory,[filename '.csv'])));
end
[nrows,ncols]=size(newsplitup);
mystr=['%s,'];
for cols=1:ncols-2
    mystr=[mystr ' ' '%s,'];
end
%mystr=[mystr '\n'];
for row=1:nrows
   fprintf(fid,[mystr '%s\n'],newsplitup{row,:});
end
fclose(fid);









