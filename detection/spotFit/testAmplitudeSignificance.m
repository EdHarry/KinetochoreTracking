function discardAmplitudeIndex = testAmplitudeSignificance(coord,covarianceMatrix,ampList,cutOff,degreeOfFreedom,debug)
%testAmplitudeSignificance does a t-test to attest if the amplitude is different from 0
%
% SYNOPSIS: discardAmplitudeIndex = testAmplitudeSignificance(coord,covarianceMatrix,cutOff,degreeOfFreedom)
%
% INPUT coord : m by n coordinate list. n is the number of dimension.
%       covarianceMatrix : covariance matrix for the amplitude only.
%       ampList : amplitude list.
%       cutOff : cut off for the T-Test
%       degreeOfFreedom : The degree of freedom of the function.
%       debug : Flag - 0,1 - if on, print additional helpful information.
% OUTPUT discardAmplitudeIndex : The index of the coordinate that failed
%        the T-Test.
%
% REMARKS
%
% SEE ALSO spotMMFit, testDistanceSignificance
%
% EXAMPLE
%

% created with MATLAB ver.: 7.14.0.739 (R2012a) on Mac OS X  Version: 10.6.8 Build: 10K549
%
% created by: Jacques Boisvert
% DATE: 11-Jun-2012
%
% Last revision $Rev: 2933 $ $Date: 2012-07-13 15:03:00 -0400 (Fri, 13 Jul 2012) $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    debug = 0;
end
nDim = size(coord,2);
amplitudeCutoff = tinv(1-cutOff,degreeOfFreedom);
nSpot = size(coord,1);
discardAmplitudeIndex = zeros(nSpot,1);
%Test every spot for zero amplitude
for cSpot = 1:nSpot
    %null-hypothesis: the amplitude = 0
    %H1: amplitude ~= 0
    diagIdx = nDim*nSpot + cSpot;
    testValue = ampList(cSpot)/sqrt(covarianceMatrix(diagIdx,cSpot));
    %if H0 accepted, we throw the spot out (we don't want spots with zero amplitude!)
    if testValue < amplitudeCutoff
        if debug
            if nDim == 2
                fprintf('Spot [%.3f %.3f] has failed amplitude t-test\n',coord(cSpot));
            elseif nDim == 3
                fprintf('Spot [%.3f %.3f %.3f] has failed amplitude t-test\n',coord(cSpot));
            end
            fprintf ('probability value : %.3f\n',testValue);
            fprintf ('cutoff : %.3f\n',amplitudeCutoff);
        end
        discardAmplitudeIndex(cSpot) = cSpot;
    end
end

discardAmplitudeIndex = find(discardAmplitudeIndex);
