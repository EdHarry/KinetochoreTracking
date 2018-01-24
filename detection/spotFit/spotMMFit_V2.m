function [fullCoordList,fullAmpList,fullBgList,covarianceMatrix,newCentroidList,stats] = spotMMFit_V2(frame,iniCoord,psfSize,varargin)
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
%
% EHarry changes (Nov 2012):
%               really major changes (big fixes):
%                                           main (this) program:
%                                                                           control of fitNPlusOne while loop passed to independent variable, not the main one in the 'para' struct, otherwise only the 1st cluster was nPlusOne fitted
%                                                                           
%                                           sub programs:
%                                                           clusterSpot:
%                                                                           to cluster in 3D with in principle different overlap factors (or cutoff) in Z compared to XY the coordinates in Z were multiplied by cutoff(1)/cutoff(3), this 'warped' the Z coordinates so that cutoff(1) could be used as before   
%                                                                           in principle this could be done for Y as well if different cutoffs were used for X compared to Y    
%                                                           
%                                                           mixtureModelFitting:
%                                                                           input solution vector to lsqnonlin has to be 1D with the inputs in the order of the columns in the jacobian, so reshapped the 'centroid' input into a 1D column vector [x;y;z;x;y;z...] before calling lsqnonlin. Subsequently the vector was reshapped back to 2D straight away in the obective function, and also the final solution 'newCoord' was reshaped to 2D after lsqnonlin.  
%
%                                               
%                                                                           For a proper sparse jacobian then the sparse pattern has to be input to the objective function (the 'JacobPattern' option is only for when finite differencing is used). The sparse pattern was calculated like before, and like in clusterSpots the Z coorinates were warped so that a single cutoff value could be used in createSpaseDistanceMatrix. The pattern was made a sparse logical matrix and was input into the objective function as a 5th variable.
%
%
%
%
%               major changes:                
%                                           main (this) program:            
%                                                                           
%                                                                           after deciding on a new centroid for nPlusOne fitting, the maskCoordinates and maskAmplitudes are re-calculated taking the new spot into account. If the new coordinates are differnt from before then the fitting is redone first with N spots then with N+1 spots (so to compare fits over the same pixels)
%                                                                           this is done to avoid the problem of adding new spots at the edge of a cluster and not having enough signal pixels to fit the new spot
%
%                                                                           after sucessfully adding a new spot to a cluster then the spot clustering is re-done to check for any clusters that have merged together because of the extra spot. If clusters have merged then the fitting is restarted with the new list of clusters
%
%
%
%
%
%               minor changes:              main (this) program:
%                                                                           position and amplitude uncertainties are also returned in the final output, the coordinates outputed go [x,y,z,dx,dy,dz] and the amplitudes go [A,dA]
%                                   
%                                                                           debug printing of spot positions that fail separation or amplitude tests will print the Z coordinate if in 3D
%                                           
%                                           testAmplitudeSignificance:      
%                                                                           debug printing of spot positions that fail amplitude tests will print the Z coordinate if in 3D
%                               
%                                           testDistanceSignificance:       debug printing of spot positions that fail separation tests will print the Z coordinate if in 3D
%
%
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
    'overlapFactor', psfSize*6 ,...% ~6 sigma
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

%% READ INTENSITIES PER CLUSTER
% will recluster if any clusters merge
reCluster = true;
while reCluster
    % initially set reCluster to false, and set to true if clusters do
    % merge
    reCluster = false;
    %Pre-allocation of exit coordinate list
    fullCoordList= [];
    fullAmpList = [];
    fullBgList = [];
    
    %Pre-allocation of coordinate uncertanties
    fullCoordUncList= [];
    fullAmpUncList = [];
    fullBgUncList = [];
    
    %debug stats
    newCentroidList = [];
    stats = [];
    for iCluster = 1:nCluster
        cCoord = clusterList{iCluster};
        
        
        %READ INTENSITIES
        if nDim == 1
            linInd = cCoord;
        elseif nDim == 2
            linInd = sub2ind(size(frame),ceil(cCoord(:,1)),ceil(cCoord(:,2)));
        elseif nDim == 3
            linInd = sub2ind(size(frame),ceil(cCoord(:,1)),ceil(cCoord(:,2)),ceil(cCoord(:,3)));
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
        warning off
        while(fitting)
            
            %Define number of fitting spots;
            nSpot = size(newCoord,1);
                        
            [newCoord,resi,newAmp,bg,jacobianMatrix] = mixtureModelFitting(newCoord,maskAmp,maskCoord,psfSize,[],size(frame));
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
        if isempty(newCoord)
            newCoord = [];
        end
        
        %% FIT N+1
        %if(para.fitNPlusOne && noRejection)
        if(para.fitNPlusOne && ~isempty(newCoord))
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
                
                %re-Prepare mask cluster
                switch nDim
                    
                    case 1
                        gMask1D = zeros(size(frame));
                        for id = 1:size(nPlusOneCoord,1)
                            gMask1D = gMask1D +  GaussMask1D(psfSize,size(frame),nPlusOneCoord(id,:),0)';
                        end
                        gMask1D = gMask1D > exp( - (para.influenceRadius^2 / 2));
                        newMaskAmp = frame(gMask1D);
                        newMaskCoord = find(gMask1D);
                        
                    case 2
                        %------------------------------2D-----------------------------------
                        %Accumulate all the gaussian sigma centered on each centroid.
                        gMask2D = zeros(size(frame));
                        for id = 1:size(nPlusOneCoord,1)
                            gMask2D = gMask2D + GaussMask2D(psfSize,size(frame),nPlusOneCoord(id,:),0,1,1);
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
                        
                        newMaskAmp = frame(gMask2D);
                        newMaskCoord = find(gMask2D);
                        [x,y] = ind2sub(size(frame),newMaskCoord);
                        newMaskCoord = [x,y];
                        
                        %-------------------------------------------------------------------
                        
                        
                        %------------------------------3D-----------------------------------
                    case 3
                        %Accumulate all the gaussian sigma centered on each centroid.
                        gMask3D = zeros(size(frame));
                        for id = 1:size(nPlusOneCoord,1)
                            gMask3D = gMask3D + GaussMask3D(psfSize,size(frame),nPlusOneCoord(id,:),0,1,1);
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
                        
                        newMaskAmp = frame(gMask3D);
                        newMaskCoord = find(gMask3D);
                        [x,y,z] = ind2sub(size(frame),newMaskCoord);
                        newMaskCoord = [x,y,z];
                        %-------------------------------------------------------------------
                end
                
                if size(newMaskCoord,1) ~= size(maskCoord,1) || ~all(newMaskCoord(:) == maskCoord(:))
                    %re-fit N with new coords if they're different
                    if size(newCoord,1) == 0
                        tmp=0;
                    end
                    [newCoord,resi,amp,bg,jacobianMatrix] = mixtureModelFitting(newCoord,newMaskAmp,newMaskCoord,psfSize,[],size(frame));
                    nFreeParam = size(jacobianMatrix,2);
                    degreeOfFreedom = size(resi,1) - nFreeParam;
                end
                
                % fit N+1
                if size(nPlusOneCoord,1) == 0
                    tmp=0;
                end
                [nPlusOneCoord,newResi,newAmp,newBg,nPlusOneJacobianMatrix] = mixtureModelFitting(nPlusOneCoord,newMaskAmp,newMaskCoord,psfSize,[],size(frame));
                
                
                %Number of free param for the n+1 fit
                newFreeParam = size(nPlusOneJacobianMatrix,2);
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
                            
                            maskAmp = newMaskAmp;
                            maskCoord = newMaskCoord;
                            
                            % test for merging clusters
                            clListNew = clusterSpot(cat(1,fullCoordList,vertcat(clusterList{iCluster+1:end}),newCoord),para.overlapFactor,0);
                            %Pre-Allocation of cell cluster
                            nClusterNew = size(unique(clListNew),1);
                            
                            if nClusterNew < nCluster % if clusers have merged then set the new cluster to the cluster list and restart
                                nCluster = nClusterNew;
                                clusterListNew = cell(nCluster,1);
                                
                                iniCoordNew = cat(1,fullCoordList,vertcat(clusterList{iCluster+1:end}),newCoord);
                                %Fill coordList for each cluster
                                for idx = 1:nCluster;
                                    clusterListNew{idx} = iniCoordNew(clListNew(:,1) == idx,:);%Could use x or y col.
                                end
                                clusterList = clusterListNew;
                                
                                % set reCluster flag, also stop fitting
                                % nPlusOne
                                fitNPlusOne = false;
                                reCluster = true;
                                
                                if para.debug
                                    fprintf('Clusters have merged. Restarting.\n');
                                end
                            end
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
        warning on
        
        % if reCluster then break for loop
        if reCluster
            break
        end
        
        %Give a background value to every coordinates.
        nSpot = size(newCoord,1);
        bg = bg * ones(nSpot,1);
        
        %Allocation of kept coordinates;
        fullCoordList = cat(1,fullCoordList,newCoord);
        fullAmpList = cat(1,fullAmpList,amp);
        fullBgList = cat(1,fullBgList,bg);
        
        if nSpot > 0
            % get uncertainties from covariance matrix
            uncertanties = sqrt(diag(covarianceMatrix));
            % take out bg uncertainty
            bgUnc = uncertanties(end);
            % assign bgUnc to each spot
            bgUnc = bgUnc * ones(nSpot,1);
            % get the amplitude uncertanties
            ampUnc = uncertanties(end-1:-1:end-nSpot);
            ampUnc = ampUnc(end:-1:1);
            % reshape the rest
            uncertanties = reshape(uncertanties(1:end-1-nSpot),nDim,nSpot)';
            
            %Allocation of kept coordinate uncertainties;
            fullCoordUncList = cat(1,fullCoordUncList,uncertanties);
            fullAmpUncList = cat(1,fullAmpUncList,ampUnc);
            fullBgUncList = cat(1,fullBgUncList,bgUnc);
        end
    end
end

% add uncertainties to the list of coords and amplitues
fullCoordList = [fullCoordList fullCoordUncList];
fullAmpList = [fullAmpList fullAmpUncList];
fullBgList = [fullBgList fullBgUncList];