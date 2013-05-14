function mat=quikrTrain(inputfasta,k)
%matrix=QuikrTrain(inputfasta,k) Returns the (sparse) k-mer sensing
%matrix from the FASTA file located at 'inputfasta' for the k-mer size k.
%This serves to retrain the Quikr method. The ouput can then be fed to
%quikrCustomTrained().
%Works for k=1:8 (typically the matrices get too large for k>9)
%The filename for the inputfasta must be the complete path.

if nargin>2
    error('too many input arguments');
end


[pathtofile,filename,ext]=fileparts(inputfasta); %get the parts of the file


outputfilename=fullfile(pathtofile, [filename sprintf('-sensingmatrixK%d.txt',k)]); %Currently this is coded to write a temporary file. In future versions, this will be all be done in RAM
%The reason for writing the file to disk first is that Matlab typically
%crashes when unix() returns as many entries as ./probabilities-by-read
%does (on the order of ~2*10^10).

kmerfilename=sprintf('%dmers.txt',k); %This contains the list of 6-mers to count. In future versions this will be computed locally instead of being read in.

unix(['probabilities-by-read ' sprintf('%d',k) ' ' inputfasta ' ' kmerfilename ' > ' outputfilename]); %obtain the k-mer counts of the inputfasta read-by-read

fid=fopen(outputfilename); %open the output file

%A=textscan(fid,'%f'); %get all the counts
%A=A{:};
A=fscanf(fid,'%f');
mat=sparse(reshape(A,4^k,length(A)/4^k)); %form into a matrix
mat=bsxfun(@rdivide,mat,sum(mat,1)); %column-normalize
fclose(fid); %close file
delete(outputfilename); %delete the file

