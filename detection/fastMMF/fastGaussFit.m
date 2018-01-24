function [coordList,ampList,bgList] = fastGaussFit( image,cands,psfSigma )
%FASTGAUSSFIT based on spotFit.m and detectSubResFeatures3D.m
%   EHarry Dec 2012

%% initialise

doMMF = 1;

% get image size
[numPixelsX,numPixelsY,numPixelsZ] = size(image);

%% first fitting

% first jointly fit all with no N+1 to remove false spots
psfSigmaFull = psfSigma([1 1 2]);
[coordList,ampList,bgList] = spotMMFit(image,cands,psfSigmaFull,'influenceRadius',1.7*psfSigmaFull,'overlapFactor',3.5*psfSigmaFull);






%% Mixture Model Fitting

% initial cluster
clusters = clusterSpots([coordList(:,1:3) ampList(:,1)],[numPixelsX,numPixelsY,numPixelsZ],psfSigma);

%set optimization options
options = optimset('Jacobian','on',...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

% default F-test cutoff
alphaR_def = 0.05;

% recluster while loop
reCluster = 1;
while reCluster
    % set initially to no reCluster
    reCluster = 0;
    % results matrix (x,y,z)
    spots = [];
    
    %go over all clusters
    iCluster = 0;
    while iCluster < length(clusters)
        
        iCluster = iCluster + 1;
        keepCluster = 1;
        
        %get initial guess of positions and amplitudes
        numMaximaT = clusters(iCluster).numMaxima;
        maximaPosT = clusters(iCluster).maximaPos(:,1:3);
        maximaAmpT = clusters(iCluster).maximaAmp;
        bgAmpT = mean(bgList(:,1));
        clusterPixels = clusters(iCluster).pixels(:,1:3);
        clusterPixels1DIdx = clusters(iCluster).pixels(:,4);
        
        %% remove superfluous maxima
        if numMaximaT > 1
            
            
            % UPDATE 120925 EHarry, a better way to find superfluous maxima
            % treating XY and Z totally separate
            
            
            % first create a 3D distance matrix
            disMat3D = createDistanceMatrix(maximaPosT,maximaPosT);
            % now separate XY and Z dis mats
            distMatXY = createDistanceMatrix(maximaPosT(:,1:2),maximaPosT(:,1:2));
            distMatZ = abs(createDistanceMatrix(maximaPosT(:,3),maximaPosT(:,3)));
            
            % now find spots above the cutoff (2*psf in XY and Z)
            distMatXY = distMatXY > 2 * psfSigma(1);
            distMatZ = distMatZ > 2 * psfSigma(2);
            
            % make the diagonals of the mats true
            for iLoop = 1:size(distMatXY,1)
                distMatXY(iLoop,iLoop) = true;
                distMatZ(iLoop,iLoop) = true;
            end
            
            % cut off mat, matrix of spots beond the cuttoff in XY OR Z
            cutOffMat = distMatXY | distMatZ;
            
            % find bad spots who don't have all trues
            badSpots = find(~all(cutOffMat, 2));
            
            % while badSpots is not empty
            while ~isempty(badSpots)
                % now find the bad spots's neibouring spots with the lowest 3D
                % distance
                minDis = Inf;
                for iSpot = badSpots'
                    for jSpot = find(~cutOffMat(iSpot,:))
                        dis = disMat3D(iSpot,jSpot);
                        if dis < minDis
                            minDis = dis;
                            candToRem = [iSpot; jSpot];
                        end
                    end
                end
                
                %determine which one of them has the smaller average distance to
                %the other maxima
                
                aveDistIJ = mean(disMat3D(candToRem,:),2);
                max2remove = candToRem(aveDistIJ==min(aveDistIJ));
                max2keep = setdiff((1:numMaximaT)',max2remove(1));
                
                %remove it from the cluster
                numMaximaT = numMaximaT - 1;
                maximaPosT = maximaPosT(max2keep,:);
                maximaAmpT = maximaAmpT(max2keep,:);
                
                % if still more than one maxima repeat
                
                if numMaximaT > 1
                    disMat3D = createDistanceMatrix(maximaPosT,maximaPosT);
                    distMatXY = createDistanceMatrix(maximaPosT(:,1:2),maximaPosT(:,1:2));
                    distMatZ = abs(createDistanceMatrix(maximaPosT(:,3),maximaPosT(:,3)));
                    distMatXY = distMatXY > 2 * psfSigma(1);
                    distMatZ = distMatZ > 2 * psfSigma(2);
                    for iLoop = 1:size(distMatXY,1)
                        distMatXY(iLoop,iLoop) = true;
                        distMatZ(iLoop,iLoop) = true;
                    end
                    cutOffMat = distMatXY | distMatZ;
                    badSpots = find(~all(cutOffMat, 2));
                    
                else
                    badSpots = [];
                end
                
            end
        end
        
        %% attempt first fit and iterative mixture-model fitting if requested
        
        %crop part of image that is relevant to this fitting
        imageC = image(clusterPixels1DIdx);
        
        firstFit = 1; %logical variable indicating first fit attempt
        fit = 1; %logical variable indicating whether to attempt to fit
        
        %fit if flag is 1
        while fit
            
            
            %clusterPixels = clusterPixels_new;
            
            
            [x0,lb,ub] = mmfInitGuessLowerUpperBounds3D(maximaPosT,maximaAmpT,...
                bgAmpT,psfSigma,clusterPixels,firstFit);
            
            %calculate number of degrees of freedom in system
            %         numDegFreeT = size(clusterPixels,1)-3*numMaximaT-1;
            numDegFreeT = size(clusterPixels,1)-4*numMaximaT-1;
            
            %determine feature positions and amplitudes and estimate background
            %intensity, using nonlinear least squares data fitting
            %also get residuals and Jacobian matrix to calculate variance-covariance
            %of estimated parameters
            [solutionT,residualsT] = fitNGaussians3D_mexCode_fitFun(options,x0, lb, ub, imageC, clusterPixels, psfSigma, [numPixelsX numPixelsY numPixelsZ], clusterPixels1DIdx);
            
            residualsT = -residualsT; %minus sign so that residuals = real image - model image
            residVarT = (sum(residualsT.^2)/numDegFreeT);
            
            %check whether addition of 1 PSF has significantly improved the fit
            if firstFit %if this is the first fit
                
                firstFit = 0; %next one won't be
                
            else %if this is not the first fit
                
                %get test statistic, which is F-distributed
                testStat = residVarT/residVar;
                
                %get p-value of test statistic
                pValue = fcdf(testStat,numDegFree,numDegFreeT);
                
                %compare p-value to alpha
                %1-sided F-test: H0: F=1, H1: F<1
                % new, if less than 4 localMaxima the nintroduce a v small
                % cutoff
                if numMaxima < 4
                    alphaR = 1e-50;
                else
                    alphaR = alphaR_def;
                end
                
                if pValue > alphaR %if p-value is larger, do not accept this fit and exit
                    fit = 0;
                else %else if fit is accepted then recluster to test for merging of clusters, also get aditional pixels that need to be fitted
                    % add fitted spots to this cluster
                    
                    clusters(iCluster).maximaPos = [solutionT(1:4:end-1) solutionT(2:4:end-1) solutionT(3:4:end-1)];
                    clusters(iCluster).maximaAmp = solutionT(4:4:end-1);
                    clusters(iCluster).numMaxima = numMaximaT;
                    % re-cluster
                    clusters_tmp = clusterSpots(clusters,[numPixelsX,numPixelsY,numPixelsZ],psfSigma);
                    % check number of clusters
                    if length(clusters_tmp) < length(clusters)% if so then keep the new cluster struct and re-start the fitting
                        clusters = clusters_tmp;
                        reCluster = 1;
                        fit = 0;
                    else % else just fit current cluster and get new pixels to fit (if any), and re-fit
                        clusters_tmp = clusterSpots(clusters(iCluster),[numPixelsX,numPixelsY,numPixelsZ],psfSigma);
                        if length(clusters_tmp) > 1
                            tmp = clusters_tmp;
                            clear clusters_tmp
                            tmp1 = catStruct(1,'tmp.pixels');
                            clusters_tmp.pixels = tmp1;
                            clear tmp1 tmp
                        end
                        if length(clusters_tmp.pixels(:,4)) ~= length(clusters(iCluster).pixels(:,4)) || ~all(clusters_tmp.pixels(:,4) == clusters(iCluster).pixels(:,4))
                            clusterPixels = clusters_tmp.pixels(:,1:3);
                            clusterPixels1DIdx = clusters_tmp.pixels(:,4);
                            imageC = image(clusterPixels1DIdx);
                            numDegFreeT = size(clusterPixels,1)-4*numMaximaT-1;
                            [solutionT,residualsT] = fitNGaussians3D_mexCode_fitFun(options,solutionT, lb, ub, imageC, clusterPixels, psfSigma, [numPixelsX numPixelsY numPixelsZ], clusterPixels1DIdx);
                            residualsT = -residualsT; %minus sign so that residuals = real image - model image
                            residVarT = (sum(residualsT.^2)/numDegFreeT);
                        end
                    end
                end
                
            end
            
            if fit %if this fit is accepted (which is the default if it's the first fit)
                
                %update variables
                numMaxima = numMaximaT;
                numDegFree = numDegFreeT;
                solution = solutionT;
                residuals = residualsT;
                residVar = residVarT;
                
                
                
                %extract estimate of background intensity and
                %remove from vectors
                bgAmp = solution(end);
                solution = solution(1:end-1);
                
                %reshape 4nx1 vector "solution"into nx3 matrices
                solution = reshape(solution,4,numMaxima);
                solution = solution';
                
                %extract feature positions and amplitudes and their uncertainties
                maximaPos = solution(:,1:3);
                maximaAmp = solution(:,4);
                
                %if attempting iterative mixture-model fitting
                if doMMF
                    
                    %add a kernel positioned at the pixel with maximum residual
                    %for fit with numMaxima + 1 Gaussians
                    numMaximaT = numMaxima + 1; %update number of maxima
                    maximaAmpT = [maximaAmp(:,1); mean(maximaAmp(:,1))]; %signal amplitude
                    indxMaxRes = find(residuals==max(residuals)); %index of pixel with maximum residuals
                    coord = clusterPixels(indxMaxRes(1),:); %position of new kernel (using indxMaxRes(1) avoids crashing if it so happens that two pixels have the exact same value which happens to be the maximum value - it has happened!!!)
                    maximaPosT = [maximaPos(:,1:3); coord]; %signal positions
                    bgAmpT = bgAmp; %background amplitude
                    
                else %if no iterative mixture model fitting
                    
                    fit = 0;
                    
                end
                
            end %(if fit)
            
        end %(while fit)
        
        % if reClustering break here
        if reCluster
            break
        end
        
        
        %% estimate of the parameters and the residuals from the fit
        if keepCluster
            
            %collect initial guesses and lower and upper bounds for fit
            %clusterPixels = clusterPixels_new;
            [x0,lb,ub] = mmfInitGuessLowerUpperBounds3D(maximaPos(:,1:3),...
                maximaAmp(:,1),bgAmp,psfSigma,clusterPixels,1);
            
            
            %calculate number of degrees of freedom in system
            numDegFree = size(clusterPixels,1)-4*numMaxima-1;
            
            %determine feature positions and amplitudes and estimate background
            %intensity, using nonlinear least squares data fitting
            %also get residuals and Jacobian matrix to calculate variance-covariance
            %of estimated parameters
            
            solution = fitNGaussians3D_mexCode_fitFun(options, x0, lb, ub, imageC, clusterPixels, psfSigma, [numPixelsX,numPixelsY,numPixelsZ], clusterPixels1DIdx);
            
            % test for merging clusters
            clusters(iCluster).maximaPos = [solution(1:4:end-1) solution(2:4:end-1) solution(3:4:end-1)];
            clusters(iCluster).maximaAmp = solution(4:4:end-1);
            clusters(iCluster).numMaxima = numMaxima;
            % re-cluster
            clusters_tmp = clusterSpots(clusters,[numPixelsX,numPixelsY,numPixelsZ],psfSigma);
            
            if length(clusters_tmp) < length(clusters)
                clusters = clusters_tmp;
                reCluster = 1;
                break
            end
            
            %extract estimate and std of background intensity and
            %remove from vectors
            solution = solution(1:end-1);
            
            %reshape 4nx1 vectors "solution" and "standDevVec" into nx4 matrices
            solution = reshape(solution,4,numMaxima);
            solution = solution';

            %extract feature positions and amplitudes and their uncertainties
            maximaPos = solution(:,1:3);
            
        end %(if keepCluster)
        
        % discard cluster if not keeping it and reset cluster count
        if ~keepCluster
            clusters(iCluster) = [];
            iCluster = iCluster - 1;
            continue
        end
        
        % add results to spots matrix
        
        spots = [spots; maximaPos]; %#ok<AGROW>
        
    end %while iCluster < length(clusters)
end % while reClusrer


%% final fit
[coordList,ampList,bgList] = spotMMFit(image,spots,psfSigmaFull,'influenceRadius',3*psfSigmaFull,'overlapFactor',1000*psfSigmaFull);


%% SUBFUNCTIONS
% function to cluster pixels and get relevent pixels for fitting
    function clusters = clusterSpots(spots,imageSize,psfSigma)
        
        % can input a cluster struct instead
        if isstruct(spots)
            spots_ = catStruct(1,'spots.maximaPos');
            spots_(:,4) = catStruct(1,'spots.maximaAmp');
            spots = spots_;
        end
        
        % get pixel locations of spots
        pos = round(spots(:,1:3));
        % get 1D locations of these
        pos1D = sub2ind(imageSize,pos(:,1),pos(:,2),pos(:,3));
        % make image of spot pixel locations
        maxPos0 = zeros(imageSize);
        maxPos0(pos1D) = 1;
        
        %generate local maximum mask
        tmpSize = round(9*psfSigma);
        tmpSize = tmpSize + (1-mod(tmpSize,2));
        template = true(tmpSize(1),tmpSize(1),tmpSize(2));
        psfRange = floor(tmpSize/2);
        
        %place local maximum masks in image
        image_ = false(imageSize(1)+2*psfRange(1),imageSize(2)+2*psfRange(1),imageSize(3)+2*psfRange(2));
        for i=1:size(spots,1)
            x0_ = pos(i,1)+psfRange(1);
            y0 = pos(i,2)+psfRange(1);
            z0 = pos(i,3)+psfRange(2);
            ymin = y0-psfRange(1);
            ymax = y0+psfRange(1);
            xmin = x0_-psfRange(1);
            xmax = x0_+psfRange(1);
            zmin = z0-psfRange(2);
            zmax = z0+psfRange(2);
            image_(xmin:xmax,ymin:ymax,zmin:zmax) = template;
        end
        image_ = image_(psfRange(1)+1:end-psfRange(1),psfRange(1)+1:end-psfRange(1),psfRange(2)+1:end-psfRange(2));
        
        %get connectivity between local maximum masks
        CC = bwconncomp(image_);
        nIsland = CC.NumObjects;
        
        % initialize clusters (SB)
        clusters(1:nIsland) = struct('numMaxima',[],'maximaPos',[],'maximaAmp',[],'pixels',[]);
        
        %find number of local maxima and their centers and amplitudes in each island
        for i=1:nIsland
            
            %get pixels making up island
            rc = CC.PixelIdxList{i};
            
            %find initial position of PSF centers in island (in pixels)
            rcCenterL = ismember(pos1D,rc(maxPos0(rc)==1));
            maximaPos_ = spots(rcCenterL,1:3);
            
            % get amplitudes of spots in cluster
            maximaAmp_ = spots(rcCenterL,4);
            
            %determine initial number of PSFs in island
            numPSF0 = size(maximaPos_,1);
            
            % get indicies to pixels in cluster
            [rcInd1,rcInd2,rcInd3] = ind2sub(imageSize,rc);
            
            %collect data in structure
            clusters(i).numMaxima = numPSF0;
            clusters(i).maximaPos = maximaPos_;
            clusters(i).maximaAmp = maximaAmp_;
            clusters(i).pixels = [rcInd1 rcInd2 rcInd3 rc];
        end
    end

end

