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
redVariation = (0:nTimepoints-2)'/(nTimepoints-2);
greenVariation = (nTimepoints-2:-1:0)'/(nTimepoints-2);
blueVariation = [(0:numTimePlotOver2-1)'/(numTimePlotOver2-1);...
    (nTimepoints-numTimePlotOver2-2:-1:0)'/(nTimepoints-numTimePlotOver2-1)];

%get the overall color per time interval
cmap = [redVariation greenVariation blueVariation];