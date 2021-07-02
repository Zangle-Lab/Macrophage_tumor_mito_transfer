function [LocList numlocs] = getloclist(froot, fstart, fext)
%function [LocList numlocs] = getloclist(froot, fstart, fext)
%function to get the list of location strings for all files within a given
%directoy
%will skip any filenames starting in the charcter '%'
%CNSIData indicates if the data format is that used at CNSI

filelist = dir([froot, fstart, '*',fext]);
fileNames = char(sort_nat({filelist.name}'));
numf = length(fileNames);

strs = cell(numf,1);

for ii = 1:numf
    S = textscan(fileNames(ii,:), '%[QPM40X]%s[frame]%s','Delimiter','_');  
    strs{ii} =char(S{2});
end

LocList = unique(strs);
LocList = LocList(~strcmp(LocList, ''));
LocList = sort_nat(LocList);

numlocs = length(LocList);