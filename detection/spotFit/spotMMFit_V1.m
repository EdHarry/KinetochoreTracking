function [fullCoordList,fullAmpList,fullBgList,covarianceMatrix,newCentroidList,stats] = spotMMFit_V1(frame,iniCoord,psfSize,varargin)
%SPOTFIT is the main function of the spot mixture model fitting
%
% SYNOPSIS: [newCoord,amp,bg,covarianceMatrix,addedCentroidList,stats] = spotMMFit(frame,iniCoord,psfSize,varargin)
%
% INPUT frame : N x M image.
%       iniCoord : N-coordinates by M-dimension.
%       psfSize  : Estimated of the spot gaussian.
%       varargin : List of variable arguments (All are optional and have
%                                              defaults)
%                 influenceRadius - Influence radius of a spot.
%                                   Default : 4 * psfSize
%                 overlapFactor - Threshold distance for clustering.
%                                 Default : 6 * psfSize
%                 fitOnlyIsolated
%                 fitNPlusOne - Flag - Try to fit N+1 spot if all the spots
%                                      are good
%                 F_test_prob - F-Test acceptance level.
%                               Default : 0.05
%                 T_test_prob - T-Test acceptance level. Distance test.
%                               Default : 0.05
%                 amplitudeCutoff - T-Test acceptance level. Amplitude
%                                   test.
%                               Default : 0.01
%                 maxIter - Number of maximum iteration for the lsqnonlin
%                           optimization
%                 debug - Flag - Print useful information.
%
% OUTPUT          coordList : final list of optimized coordinates
%                 amp       : list of amplitude for each spots
%                 bg        : list of background intensities for each spots
%          covarianceMatrix : The Covariance Matrix used for
%                             intensity and distance statistical test.
%                stats      : Vector containing :
%                             chi1,nFreeParameter,chi2,nFreeParameter,
%                             size of the mask
%                             coord
%
% REMARKS
%
% SEE ALSO clusterSpot , mixtureModelFitting
%
% EXAMPLE
%

% created with MATLAB ver.: 7.14.0.739 (R2012a) on Mac OS X  Version: 10.6.8 Build: 10K549
%
% created by: Jacques Boisvert
% DATE: 24-May-2012
%
% Last revision $Rev: 3153 $ $Date: 2012-11-09 13:45:03 -0500 (Fri, 09 Nov 2012) $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------
%% VALIDATION OF MENDATORY INPUTS
%Define number of dimension -- 2 || 3
if size(frame,3) > 1
    nDim = 3;
elseif size(frame,2) > 1 && size(frame,1) > 1
    nDim = 2;
else
    nDim = 1;
end

if nDim == 1
    if size(frame,2) > 1
        frame = frame';
    end
end
if isempty(iniCoord)
    return;
end


%% DEFAULT ASSIGNMENT OF VARIABLE ARGUMENTS

para = struct('influenceRadius' ,psfSize*4,... ~ 4sigma
    'overlapFactor', psfSize*9 ,...% ~6 sigma
    'fitOnlyIsolated' , 0 , ...
    'fitNPlusOne',0,...
    'F_test_prob',[],...
    'T_test_prob',0.05,...
    'amplitudeCutoff',0.01,...
    'maxIter',[],...
    'debug',0);

%threshold = [0.347865369258613,0.237114023885974,0.107682402972059,0.028671224650176,0.004769039410945,0.003147930698649,1.575576248666999e-04,7.292868055914780e-05,1e-06];
coefficientFit = [1.0788,-0.4655,-0.0232];

fn = fieldnames(para);
%I loop through the variables names given in Parameters.
for i = 1:2:length(varargin)%(variables names are odd value of varargin
    if validatestring(varargin{i},fn)
        para.(varargin{i}) = varargin{i+1};
    end%If i find the value i replace the default value.
end

%% SPOT CLUSTERING

clList = clusterSpot(iniCoord,para.overlapFactor,0);
%Pre-Allocation of cell cluster
nCluster = size(unique(clList),1);
clusterList = cell(nCluster,1);

%Fill coordList for each cluster
for idx = 1:nCluster;
    clusterList{idx} = iniCoord(clList(:,1) == idx,:);%Could use x or y col.
end

%Pre-allocatio of exit coordinate list
fullCoordList= [];
fullAmpList = [];
fullBgList = [];
newCentroidList = [];
stats = [];
%% READ INTENSITIES PER CLUSTER
for iCluster = 1:nCluster
    cCoord = clusterList{iCluster};
    
    
    %READ INTENSITIES
    if nDim == 1
        linInd = cCoord;
    elseif nDim == 2
        linInd = sub2ind(size(frame),round(cCoord(:,1)),round(cCoord(:,2)));
    elseif nDim == 3
        linInd = sub2ind(size(frame),round(cCoord(:,1)),round(cCoord(:,2)),round(cCoord(:,3)));
    end
    newAmp = frame(linInd);
    
    
    %% KERNEL INITIATION
    %Prepare mask cluster
    switch nDim
        
        case 1
            gMask1D = zeros(size(frame));
            for id = 1:size(cCoord,1)
                gMask1D = gMask1D +  GaussMask1D(psfSize,size(frame),cCoord(id,:),0)';
            end
            gMask1D = gMask1D > exp( - (para.influenceRadius^2 / 2));
            maskAmp = frame(gMask1D);
            maskCoord = find(gMask1D);
            
        case 2
            %------------------------------2D-----------------------------------
            %Accumulate all the gaussian sigma centered on each centroid.
            gMask2D = zeros(size(frame));
            for id = 1:size(cCoord,1)
                gMask2D = gMask2D + GaussMask2D(psfSize,size(frame),cCoord(id,:),0,1,1);
            end
            if size(para.influenceRadius,2) == 1
                para.influenceRadius = [para.influenceRadius,para.influenceRadius];
            end
            %Threshold of the influence radius is 4*sigma;
            %Gaussian =exp(-x^2/2*sigma^2) if x = 4*sigma -> exp(-8);
            % take average of influenceFactors in case we make it not a scalar
            % some day
            influenceFactor = (mean(para.influenceRadius./psfSize)).^2/2;
            gMask2D = gMask2D > exp( - influenceFactor );
            
            maskAmp = frame(gMask2D);
            maskCoord = find(gMask2D);
            [x,y] = ind2sub(size(frame),maskCoord);
            maskCoord = [x,y];
            
            %-------------------------------------------------------------------
            
            
            %------------------------------3D-----------------------------------
        case 3
            %Accumulate all the gaussian sigma centered on each centroid.
            gMask3D = zeros(size(frame));
            for id = 1:size(cCoord,1)
                gMask3D = gMask3D + GaussMask3D(psfSize,size(frame),cCoord(id,:),0,1,1);
            end
            %Threshold of the influence radius of nSigma is n^2/2;
            %Exemple for 4*sigma : Gaussian =exp(-x^2/2*sigma^2) if x = 4*sigma -> exp(-8);
            if size(para.influenceRadius,2) == 1
                para.influenceRadius = [para.influenceRadius,para.influenceRadius,para.influenceRadius];
            end
            % take average of influenceFactors in case we make it not a scalar
            % some day
            influenceFactor = (mean(para.influenceRadius./psfSize)).^2/2;
            gMask3D = gMask3D > exp( - influenceFactor );
            
            maskAmp = frame(gMask3D);
            maskCoord = find(gMask3D);
            [x,y,z] = ind2sub(size(frame),maskCoord);
            maskCoord = [x,y,z];
            %-------------------------------------------------------------------
    end
    
    
    fitting = 1;%Bool - Fitting until all the spot are done.
    noRejection = 1; %Bool - if no spot was discarded on first run possible to try n+1 fitting
    stable = 0; %Bool - use to know if fitting has reach a stable point.
    newCoord = cCoord;
    
    %% FIT SPOTS
    while(fitting)
        
        %Define number of fitting spots;
        nSpot = size(newCoord,1);
        
        [newCoord,resi,newAmp,bg,jacobianMatrix] = mixtureModelFitting(newCoord,maskAmp,maskCoord,psfSize);
        %% CALC COVARIANCE MATRIX FROM JACOBIANMATRIX
        %Covariance is the inverse of : (jacobian^t * jacobian).
        nFreeParam = size(jacobianMatrix,2);
        degreeOfFreedom = size(resi,1) - nFreeParam;
        covarianceMatrix = inv(jacobianMatrix' * jacobianMatrix) .* (sum(resi.^2) / degreeOfFreedom);
        %% CALC STD OF BACKGROUND AND AMPLITUDE
        %% CHECK COORD OUT OF BOUND.
        
        %% TEST DISTANCE
        [discardDistanceIdx,conflictIdx] = testDistanceSignificance(newCoord,covarianceMatrix,para.T_test_prob,degreeOfFreedom,para.debug);
        
        %% TEST AMPLITUDE
        %Amplitude test
        %Get good list of spot for amplitude test.
        spotIdx = (1:nSpot)';
        if ~(isempty(conflictIdx))
            %goodIdx are the inverse of the conflictIndexes.
            goodIdx = spotIdx(~logical(sum(bsxfun(@eq,spotIdx,conflictIdx'),2)));
        else
            %if conflictIdx = [], above line will give [] when it should
            %give
            goodIdx = spotIdx;
        end
        
        %Pre-Allocation
        discardAmplitudeIdx = zeros(nSpot,1);
        %Run amplitude test
        discardFlag = testAmplitudeSignificance(newCoord(goodIdx,:),covarianceMatrix(:,nSpot*nDim+1:end-1),newAmp,para.amplitudeCutoff,degreeOfFreedom);
        if ~(isempty(discardFlag))
            discardAmplitudeIdx(goodIdx(discardFlag)) = goodIdx(discardFlag);
        end
        %Discard all the zeros from preallocation
        discardAmplitudeIdx = find(discardAmplitudeIdx);
        
        
        if ~(isempty(discardDistanceIdx)) || ~(isempty(discardAmplitudeIdx))
            noRejection = 0;
        else%All spot passed!
            stable = 1;
            if para.debug
                fprintf('No spot was rejected\n');
            end
        end
        %% IF GOLD <EJECT>
        if noRejection
            fitting = 0;
        elseif stable
            fitting = 0;
            %% ELSE REMOVE BAD SPOT
        else
            %If we have some distance rejection and amplitude rejection.
            %Start with the distance.
            if ~isempty(discardDistanceIdx) && ~isempty(discardAmplitudeIdx)
                if para.debug
                    fprintf('Distance and Amplitude test breach\n');
                    fprintf('Spot %d has not passed the distance t-test and will be discarded\n',discardDistanceIdx);
                end
                newCoord(discardDistanceIdx,:) = [];
                newAmp(discardDistanceIdx) = [];
            else
                if para.debug
                    if ~isempty(discardDistanceIdx)
                        fprintf('Spot %d has not passed the distance t-test and will be discarded\n',discardDistanceIdx);
                    else
                        fprintf('Spot %d has not passed amplitude t-test and will be discarded\n',discardAmplitudeIdx);
                    end
                    
                end
                newCoord(discardDistanceIdx,:) = [];
                newCoord(discardAmplitudeIdx,:) = [];
                newAmp(discardDistanceIdx) = [];
                newAmp(discardAmplitudeIdx) = [];
            end
            if size(newCoord,1) == 0
                fitting = 0;
            end
        end
    end
    %This prevent a bug with the concatenation of an Empty matrix : 1-by-0
    %with an array containing something.
    if isempty(newAmp)
        amp = [];
    else
        amp = newAmp;
    end
    
    %% FIT N+1
    if(para.fitNPlusOne && noRejection)
        if para.debug
            fprintf('Entering N+1 fit\n');
        end
        fitNPlusOne = true;
        while(fitNPlusOne)
            %-----Easy way to get a new set of coordinates-----
            %  take max of residual as a new possible centroid
            [~,maxIdx] = max(resi);
            newCentroid = maskCoord(maxIdx,:);
            nPlusOneCoord = cat(1,newCentroid,newCoord);
            
            if para.debug
                newCentroidList = cat(1,newCentroidList,newCentroid);
            else
                newCentroidList = [];
            end
            
            [nPlusOneCoord,newResi,newAmp,newBg,nPlusOneJacobianMatrix] = mixtureModelFitting(nPlusOneCoord,maskAmp,maskCoord,psfSize);
            %Number of free param for the n+1 fit
            newFreeParam = nFreeParam + nDim + 1;
            newDegreeOfFreedom = size(newResi,1) - newFreeParam;
            chi1 = sum(resi(:).^2)/degreeOfFreedom;
            chi2 = sum(newResi(:).^2)/newDegreeOfFreedom;
            
            
            %%%%%%% PSEUDO WORK%%%%%%%%%%
            %% F-TEST
            fValue = chi1/chi2;
            probability = fcdf(fValue,degreeOfFreedom,newDegreeOfFreedom);
            
            %The f-test threhold was obtained by taking the 5% of empirical
            %thresholds.
            %f-test value with different snr.
            if isempty(para.F_test_prob)
                %SNR is the maximum intensity divided by the square root of
                %the standard deviation of the intensity.
                
                if nDim == 2
                    snr = 0.6 * ( max(frame(sub2ind(size(frame),maskCoord(:,1),maskCoord(:,2)))) / std(frame(sub2ind(size(frame),maskCoord(:,1),maskCoord(:,2)))));
                elseif nDim == 3
                    snr = 0.6 * ( max(frame(sub2ind(size(frame),maskCoord(:,1),maskCoord(:,2),maskCoord(:,3)))) / std(frame(sub2ind(size(frame),maskCoord(:,1),maskCoord(:,2),maskCoord(:,3)))));
                end
                
                F_test_prob = coefficientFit(1) * exp(coefficientFit(2) * snr) + coefficientFit(3);
                if(F_test_prob > 0)
                    F_test_prob = 1 - F_test_prob;
                else
                    F_test_prob = 0.0000001;
                end
            end
            %%%%%%% PSEUDO WORK %%%%%%%%%%
            if para.debug
                st(1,1) = chi1;
                st(1,2) = nFreeParam;
                st(1,3) = chi2;
                st(1,4) = newFreeParam;
                st(1,5) = size(maskCoord,1);
                stats = cat(1,stats,st);
            end
            
            if(probability > F_test_prob)
                %Covariance is the inverse of : (jacobian^t * jacobian).
                nPlusOneCovarianceMatrix = inv(nPlusOneJacobianMatrix' * nPlusOneJacobianMatrix) .* (sum(newResi.^2) / newDegreeOfFreedom);
                
                %Test significance with amplitude and distance
                %% TEST DISTANCE N + 1
                [discardDistanceIdx,conflictIdx] = testDistanceSignificance(nPlusOneCoord,nPlusOneCovarianceMatrix,para.T_test_prob,newDegreeOfFreedom,para.debug);
                if ~(isempty(discardDistanceIdx))
                    if para.debug
                        fprintf('fit N+1 distance test failed\n');
                    end
                    fitNPlusOne = false;
                else %Distance Test Success
                    %% TEST AMPLITUDE N + 1
                    
                    nPlusOneSpot = size(nPlusOneCoord,1);
                    %Get the spots indexes that weren't problematic in term
                    %of distance. This list of spot will be use
                    %for amplitude test.
                    spotIdx = (1:nPlusOneSpot)';
                    if ~(isempty(conflictIdx))
                        %goodIdx are the inverse of the conflictIndexes.
                        goodIdx = spotIdx(~logical(sum(bsxfun(@eq,spotIdx,conflictIdx'),2)));
                    else
                        %if conflictIdx = [], above line will give [] when it should
                        %give
                        goodIdx = spotIdx;
                    end
                    %preAllocation
                    discardAmplitudeIdx = zeros(nPlusOneSpot,1);
                    %run amplitude test
                    discardFlag = testAmplitudeSignificance(nPlusOneCoord(goodIdx,:),nPlusOneCovarianceMatrix(:,nPlusOneSpot*nDim+1:end-1),newAmp,para.amplitudeCutoff,newDegreeOfFreedom);
                    
                    if ~(isempty(discardFlag))
                        discardAmplitudeIdx(goodIdx(discardFlag)) = goodIdx(discardFlag);
                    end
                    
                    %Discard all the zeros from preallocation
                    discardAmplitudeIdx = find(discardAmplitudeIdx,1);
                    
                    if ~(isempty(discardAmplitudeIdx))
                        if para.debug
                            fprintf('Fit N+1 amplitude test failed\n');
                        end
                        fitNPlusOne = false;
                    end
                    
                    if fitNPlusOne
                        %F-Test,distance and amplitude test have all passed.
                        %Store new coords,amp,bg,degreeOfFreedom,nFreeParam.
                        if para.debug
                            if nDim == 2
                                fprintf('Fit-N+1 pass, adding : %.2f , %.2f\n', nPlusOneCoord(1,1),nPlusOneCoord(1,2));
                            elseif nDim == 3
                                fprintf('Fit-N+1 pass, adding : %.2f , %.2f , %.2f\n', nPlusOneCoord(1,1),nPlusOneCoord(1,2),nPlusOneCoord(1,3));
                            end
                        end
                        newCoord = nPlusOneCoord;
                        amp = newAmp;
                        bg = newBg;
                        resi = newResi;
                        nFreeParam = newFreeParam;
                        degreeOfFreedom = newDegreeOfFreedom;
                        
                        covarianceMatrix = nPlusOneCovarianceMatrix;
                    end
                    
                end
            else%F-TEST fail
                if para.debug
                    fprintf('Fit N+1 F-Test failed at %.2f \n',probability);
                end
                fitNPlusOne = false;
            end
        end
    end
    %Give a background value to every coordinates.
    nSpot = size(newCoord,1);
    tmpBg = NaN(nSpot,1);
    for row = 1:nSpot
        tmpBg(row) = bg;
    end
    bg = tmpBg;
    %Allocation of kept coordinates;
    fullCoordList = cat(1,fullCoordList,newCoord);
    fullAmpList = cat(1,fullAmpList,amp);
    fullBgList = cat(1,fullBgList,bg);
end
