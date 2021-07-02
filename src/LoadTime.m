function [time] = LoadTime(fname)
%function [time] = LoadTime(fname)
%function to load the time associated with a given filename
%this time will be a datestring, as returned by the dir function
%fname is the .mat file name associated with the .tif file of the original
%data, so time is the modification time of the original tif file.

D = dir([fname(1:end-4), '.tif']);

time = D.date;

end

