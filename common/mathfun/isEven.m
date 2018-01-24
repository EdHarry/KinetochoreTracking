function out = isEven(in)
%ISEVEN checks whether a number is even
%
% SYNOPSIS out = isEven(in)
%
% INPUT    in :  input (array) of numbers to be tested. IsEven only works
%                   properly for integers. 
% OUTPUT   out:  array of size(in) with 
%                   1 for even integers and zero
%                   0 for odd integers
%                 NaN for non-integers
%
% c: jonas 5/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% out = remainder of input + 1 divided by two
out = mod(in + 1,2);

% set NaNs for not 0 or 1
out(find(~(out==1 | out==0))) = NaN;

% make logical
out = logical(out);