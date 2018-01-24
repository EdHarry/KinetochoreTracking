function [tiledArray]=tileSquaresWithIndex(im,squareLength)
%TILESQUARESWITHINDEX gives tiledArray where each entry contains tile index
%
% DESCRIPTION: tileSquaresWithIndex divides an array of size im into an
% integer number of squares with side length squareLength. If size(im) does
% not divide evenly, the output array is smaller than im. The output,
% tileArray, is simply a grid of squares where each entry in the matrix is
% the index of that particular square.
%
% SYNOPSIS: [tiledArray]=tileSquaresWithIndex(im,squareLength)
%
% INPUT: im           : image or matrix (only needed for dimensions)
%        squareLength : dimension of tiles to be created
%
% OUTPUT: tiledArray  : tiled array where the array values are the index of
%                       the particular tile
%
% e.g.) [tiledArray]=tileSquaresWithIndex(ones(10),2) gives:
% 
%      1     1     6     6    11    11    16    16    21    21
%      1     1     6     6    11    11    16    16    21    21
%      2     2     7     7    12    12    17    17    22    22
%      2     2     7     7    12    12    17    17    22    22
%      3     3     8     8    13    13    18    18    23    23
%      3     3     8     8    13    13    18    18    23    23
%      4     4     9     9    14    14    19    19    24    24
%      4     4     9     9    14    14    19    19    24    24
%      5     5    10    10    15    15    20    20    25    25
%      5     5    10    10    15    15    20    20    25    25
%
% MATLAB VERSION (originally written on): 7.2.0.232 (R2006a) Windows_NT
%
% USERNAME: kathomps
% DATE: 13-Jun-2007
%

nBr=floor(size(im,1)/squareLength); % number of boxes along rows
nBc=floor(size(im,2)/squareLength); % number of boxes along columns
Ms=reshape(1:nBr*nBc,nBr,[]);

tiledArray = round(squareLength^2*imResample(Ms, [squareLength squareLength], [1 1]));



