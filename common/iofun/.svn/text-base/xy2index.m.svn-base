function [indexList] = xy2index(x,y,L,W,unit)
%
% DESCRIPTION: Converts xy-coordinates to the array index for a matrix
%
% SYNOPSIS: [indexList] = xy2index(x,y,L,W,unit)
%
% INPUT: 
%        x and y    : vectors containing x and y coordinates
%        L and W    : length and width of matrix 
%        unit (opt) : default is 1
%                     if an image is subsampled, unit corresponds to the
%                     size of the original pixel in the new units. eg) if 1
%                     pixel is divided into 20 smaller boxes and you want
%                     to know which original pixel each of the boxes came
%                     from, you would use unit=20
%
% OUTPUT: 
%        indexList  : array indices for the matrix
%
% MATLAB VERSION (originally written on): 7.0.1.24704 (R14) Service Pack 1 Windows_NT
%
% USERNAME: kathomps
% DATE: 12-Apr-2006
%

if nargin<5
    unit=1;
end

col=ceil(x/unit);
row=ceil(y/unit);
indexList=L*(col-1)+row;
indexList(row<1 | row>L | col<1 | col>W) = nan;