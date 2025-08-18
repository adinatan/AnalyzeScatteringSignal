function [xyz, AtomNames] = AnXYZreader(filename)
% Fast reader for XYZ coordinates and atomic names from text files
% formatted as follows:
%
% NumOfAtoms
%   
% Atom1    x y z 
% Atom2    x y z 
% .....    .....
% .....    .....
% NumOfAtoms
%  
% Atom1    x y z 
% Atom2    x y z 
% .....    .....
% .....    .....
%
% with the assumption that the xyz file content can be contained in memory
% Example:

% filename= 'D:\Dropbox\Matlab\LCLS-LV96\qmmm_15shells\ptp_h2o_d3blyp2023\es\ptp_h2o_00.xyz'
% [xyz, AtomNames] = AnXYZreader(filename)

%
%   Version 1.1
%   Adi Natan (natan@stanford.edu)



% Read raw file
fid = fopen(filename, 'r');
rawData = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
rawData=rawData{1};
fclose(fid);
% Parse file:
N = length(rawData);
NumOfAtoms =  str2double(rawData{1});   
frames=N/(NumOfAtoms+2); %  skip 2 cells that contain # of atoms and an empty cell
AtomNames = textscan( strjoin(rawData(3:NumOfAtoms+2), '\n'), '%s %*f %*f %*f');
AtomNames=AtomNames{1};
toKeep=find(repmat([0 ; 0; ones(NumOfAtoms,1)],frames,1)); % the data to keep
% Concatenate all strings into a single large string separated by new lines
bigString = strjoin(rawData(toKeep(:)), '\n');
% Parse the large string
xyzArray = textscan(bigString, '%*s %f %f %f');
xyz=permute(reshape(cell2mat(xyzArray), [NumOfAtoms, frames, 3]),[1 3 2]);
end
