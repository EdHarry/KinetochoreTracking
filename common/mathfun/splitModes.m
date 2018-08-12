function cutValue = splitModes(data,loHi,jumpMax,verbose)
%SPLITMODES cut histogram where there is no or little data
%
% SYNOPSIS cutValue = splitModes(data,loHi,jumpMax)
%
% INPUT    data : data to split
%          loHi : (opt) number of minima that should be checked to the left
%                 and right of the initial guess for the cutoff. The lowest
%                 minimum right of the highest peak will be chosen.
%                 Default: [1,2] - 1 to the left, 2 to the right
%          jumpMax : (opt) if true, minimum is searched to the left of the
%                 highest maximum, too. Default: 0
%          verbose : (opt) if 1, output will be plotted. Default: 0
%
% c: 07/2008 Jonas Dorn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(loHi)
    loHi = [1,2];
end
if nargin < 3 || isempty(jumpMax)
    jumpMax = false;
end
if nargin < 4 || isempty(verbose)
    verbose = false;
end

% first guess via cutFirstHistMode - put as 1 to see histograms
[cutIdx, cutVal,sp] = cutFirstHistMode(data,0);

% now check the local minima in the vicinity of the cutoff
spder = fnder(sp);
zeroList = fnzeros(spder);
zeroList = zeroList(1,:);
% evaluate
zeroVals = fnval(sp,zeroList);

% look in zeroList. Find one value before cutVal, three after. Go into
% zeroVals and find lowest minimum
[dummy,closestIdx] = min(abs(zeroList - cutVal));

% check only the minimas that are close by; two to the right and
% one to the left (don't forget that between two minima there will
% always be a maximum!)
indexList = (closestIdx-loHi(1)*2):(closestIdx + loHi(2)*2);
indexList(indexList < 1 | indexList > length(zeroVals)) = [];

% also, never allow going left of the highest maximum!
[maxVal,maxValIdx] = max(zeroVals);
if ~jumpMax
    indexList(indexList < maxValIdx) = [];
end

% find lowest
[dummy, cutIdx] = min(zeroVals(indexList));
% and determine break value
cutValue = zeroList(indexList(cutIdx));

% if verbose: plot.
if verbose
    figure
    ah = gca;
    histogram(ah,data,1,0);
    hold on
    plot(ah,[cutVal;cutVal],[0,maxVal],':r')
    plot(ah,[cutValue;cutValue],[0,maxVal],'r')
end