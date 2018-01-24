function [discardIndex,conflictIndex] = testDistanceSignificance(coord,covarianceMatrix,cutoff,degreeOfFreedom,debug)
%TESTDISTANCESIGNIFICANCE test the distance significance of spots. H0 : distance is equal to 0 vs H1 distance is different of 0
%
% SYNOPSIS: [discardIndex,conflictIndex] testDistanceSignificance(coord,covarianceMatrix,cutoff,degreeOfFreedom)
%
% INPUT coord : nSpot-by-nDim list of coordinates
%       covarianceMatrix : full covariance matrix of the distance
%       cutoff : Significance level of the t-test
%       degreeOfFreedom : degree of freedom of the fitS
%       debug : flag - 0,1 - if on, print additional information
% OUTPUT discardIndex : The index that should be removed
%        conflictIndex : All the coordinate index that failed the distance test
% REMARKS
%
% SEE ALSO spotMMFit, testDistanceAndAplitudes
%
% EXAMPLE
%

% created with MATLAB ver.: 7.14.0.739 (R2012a) on Mac OS X  Version: 10.6.8 Build: 10K549
%
% created by: Jacques Boisvert
% DATE: 07-Jun-2012
%
% Last revision $Rev: 2933 $ $Date: 2012-07-13 15:03:00 -0400 (Fri, 13 Jul 2012) $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    debug = 0;
end

tTestCutoff = tinv((1-cutoff),degreeOfFreedom);
nSpots = size(coord,1);
nDim = size(coord,2);

%Pre-allocation
%Possible pair is n(n-1)/2
indistinguishable = zeros(2,nSpots*(nSpots-1)/2);
%Index for indistinguishable
indistinguishableIdx = 1;
for cSpot = 1:nSpots-1
    for otherSpot = cSpot+1:nSpots
        
        %covariates = blkdiag(covarianceMatrix( (cSpot-1)*nDim + 1:(cSpot-1)*nDim + (nDim),(cSpot-1)*nDim+1:(cSpot-1)*nDim+(nDim)),...
        %    covarianceMatrix( (otherSpot-1)*nDim+1:(otherSpot-1)*nDim+nDim, (otherSpot-1)*nDim+1:(otherSpot-1)*nDim+nDim));
        
        covariates = [covarianceMatrix( (cSpot-1)*nDim + 1:(cSpot-1)*nDim + (nDim),(cSpot-1)*nDim+1:(cSpot-1)*nDim+(nDim)) , covarianceMatrix( (cSpot-1)*nDim + 1:(cSpot-1)*nDim + (nDim),(otherSpot-1)*nDim+1:(otherSpot-1)*nDim+(nDim)) ; covarianceMatrix( (otherSpot-1)*nDim+1:(otherSpot-1)*nDim+nDim, (cSpot-1)*nDim+1:(cSpot-1)*nDim+nDim) , covarianceMatrix( (otherSpot-1)*nDim+1:(otherSpot-1)*nDim+nDim, (otherSpot-1)*nDim+1:(otherSpot-1)*nDim+nDim)];
        
        
        %do the standard gaussian error propagation
        difference=(coord(cSpot,:)-coord(otherSpot,:));
        distance = sqrt(sum(difference.^2));
        H=1/distance*[difference -difference];
        Qdd=H*covariates*H';
        
        sigD =sqrt(Qdd);
        testValue = distance/sigD;
        
        %test: H0: zero distance
        %H1: nonzero distance (positive->one-sided test)
        %keep only coords that pass as nonzero distance
        
        if testValue<tTestCutoff %one-sided
            if debug
                if nDim == 2
                    fprintf('Spot [%.2f %.2f] in conflict with Spot [%.2f %0.2f]\n',coord(cSpot,1),coord(cSpot,2),coord(otherSpot,1),coord(otherSpot,2));
                elseif nDim == 3
                    fprintf('Spot [%.2f %.2f %.2f] in conflict with Spot [%.2f %0.2f %0.2f]\n',coord(cSpot,1),coord(cSpot,2),coord(cSpot,3),coord(otherSpot,1),coord(otherSpot,2),coord(otherSpot,3));
                end
                fprintf('probability value : %.3f\n t-test cutoff : %.3f\n',testValue,tTestCutoff);
            end
            indistinguishable(1,indistinguishableIdx) = cSpot;
            indistinguishable(2,indistinguishableIdx) = otherSpot;
        end
        indistinguishableIdx = indistinguishableIdx+1;
    end
end

%Taking out the extra spaces define by zeros
gIdx = find(indistinguishable);
indistinguishable = reshape(indistinguishable(gIdx),2,size(gIdx,1)/2);

if~(isempty(indistinguishable))
    
    %count which spot appears the most
    [uniqueEntries,numberOfOccurences] = countEntries(indistinguishable);
    spotCount = [uniqueEntries,numberOfOccurences];
    spotCount = -sortrows(-spotCount,2); %sort in descending order
    conflictIndex = spotCount(:,1);
    %spotCount(1,1) is the spot with the most number of conflict
    %spotCount(1,2) is the maximum number of conflict.
    nEqual = size(find(spotCount(:,2) == spotCount(1,2)),1);
    %If there is no spot with more conflict than others :
    if nEqual > 1
        covTrace = zeros(nEqual,1);
        for idx = 1:nEqual
            cSpot = spotCount(idx,1);
            covTrace(idx) = trace(covarianceMatrix( (cSpot-1)*nDim + 1:(cSpot-1)*nDim + (nDim),(cSpot-1)*nDim+1:(cSpot-1)*nDim+(nDim)));
        end
        [~,maxIdx] = max(covTrace);
        discardIndex = maxIdx;
        discardIndex = conflictIndex(discardIndex);
        
    else
        %Discard worst conflicted spots.
        discardIndex = spotCount(1,1);
    end
else
    discardIndex = [];
    conflictIndex = [];
end