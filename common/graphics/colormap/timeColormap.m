function cmap = timeColormap(nTimepoints)
%TIMECOLORMAP is Khuloud's time color map
%
% SYNOPSIS: cmap = timeColormap(nTimepoints)
%
% INPUT nTimepoints: (opt) length of colormap. Default: 64
%
% OUTPUT cmap : nTimepoints-by-3 colormap that changes from green to blue to red
%
% REMARKS This code has been lifted out of Khuloud's plotTracks2D
%         In plotTracks2D, the colormap returns one fewer entry than
%         timepoints, because tracks are drawn between timepoints. This has
%         been changed.
%
% created with MATLAB ver.: 7.6.0.324 (R2008a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 16-May-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input
if nargin < 1 || isempty(nTimepoints)
    nTimepoints = 64;
end

% make colormap

%get the fraction of each color in each time interval to be plotted
numTimePlotOver2 = ceil((nTimepoints-1)/2); %needed to change blue color over time
redVariation = (0:nTimepoints-1)'/(nTimepoints-1);
greenVariation = (nTimepoints-1:-1:0)'/(nTimepoints-1);
blueVariation = [(0:numTimePlotOver2-1)'/(numTimePlotOver2);...
    (nTimepoints-numTimePlotOver2-1:-1:0)'/(nTimepoints-numTimePlotOver2)];

%get the overall color per time interval
cmap = [redVariation greenVariation blueVariation];