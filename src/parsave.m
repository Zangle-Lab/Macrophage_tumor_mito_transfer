function [ ] = parsave( fname, var1, varname, append )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

eval([varname '=var1;']);

if append
    save(fname, '-append', varname);
else
    save(fname, varname);
end

end

