function [detectedFeatures,imageN3,errFlag] = ...
    detectSubResFeatures3D_mlSparse_reCluster(image,cands,psfSigma,testAlpha,visual,...
    doMMF,bitDepth,saveResults,bgNoiseSigma)

% Edit of dectectSubResfeatures2D_V2 to fit in 3D
% EHarry March 2012

%% ORIGINAL HEADER
% % %DETECTSUBRESFEATURES2D_V2 determines the positions and intensity amplitudes of sub-resolution features using mixture model fitting
% % %
% % %SYNOPSIS [detectedFeatures,clustersMMF,imageN3,errFlag] = ...
% % %    detectSubResFeatures2D_V2(image,cands,psfSigma,testAlpha,visual,...
% % %    doMMF,bitDepth,saveResults,bgNoiseSigma)
% % %
% % %INPUT  image      : Image being analyzed.
% % %       cands      : Cands structure as output from fsmCenter.
% % %       psfSigma   : Standard deviation of point spread function (in pixels).
% % %       testAlpha  : Alpha-values for statistical tests. Structure with fields:
% % %             .alphaR: For the residuals test, comparing N+1-kernel fit to
% % %                      N-kernal fit. Optional. Default: 0.05.
% % %             .alphaA: For amplitude test. Optional. Default: 0.05.
% % %             .alphaD: For distance test. Optional. Default: 0.05.
% % %             .alphaF: Final residuals test, comparing residuals from final
% % %                      fit to estimated background noise.
% % %                      Optional. Default: 0.05.
% % %       visual     : 1 if user wants to view results; 0 otherwise.
% % %                    Optional. Default: 0.
% % %       doMMF      : 1 if user wants to do mixture-model fitting, 0
% % %                    otherwise. Optional. Default: 1.
% % %       bitDepth   : Camera bit depth. Optional. Default: 14.
% % %       saveResults: 1 if results are to be saved (in file 'detectedFeatures.mat'),
% % %                    0 otherwise. Optional. Default: 0.
% % %       bgNoiseSigma:Standard deviation of background noise. Optional. If
% % %                    not input, the code will estimate it from the image.
% % %
% % %       All optional variables can be entered as [] to use default values.
% % %
% % %OUTPUT detectedFeatures: Structure with fields:
% % %             .xCoord    : Image coordinate system x-coordinate of detected
% % %                          features [x dx] (in pixels).
% % %             .yCoord    : Image coorsinate system y-coordinate of detected
% % %                          features [y dy] (in pixels).
% % %             .amp       : Amplitudes of PSFs fitting detected features [a da].
% % %       clustersMMF: Array of clusters of sub-resolution features.
% % %                    Structure with fields:
% % %             .position  : Position of each feature in image coordinate
% % %                          system (in pixels): [x y dx dy] = [y x dy dx]
% % %                          in matrix coordinate system.
% % %             .amplitude : Intensity of each feature [A dA].
% % %             .bgAmp     : Background intensity [Abg dAbg].
% % %             .varCovMat : Variance-covariance matrix of estimated parameters.
% % %       imageN3    : Image with labeled features. Blue: those from cands;
% % %                    Red: those from mixture-model fitting; Magenta: those
% % %                    from MMF which coincide with those from cands.
% % %                    Will be output only if visual = 1.
% % %       errFlag    : 0 if function executes normally, 1 otherwise.
% % %
% % %Khuloud Jaqaman, August 2005

%% Output

detectedFeatures = [];
clustersMMF = [];
imageN3 = [];
errFlag = 0;

%% Input

%check whether correct number of input arguments was used
if nargin < 3
    disp('--detectSubResFeatures2D_V2: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check whether statistical test alpha values were inputted
if nargin < 4 || isempty(testAlpha) %if not, assign defaults
    
    testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0.05);
    
else %if some were, check their values and assign default for the rest
    
    if ~isfield(testAlpha,'alphaR')
        testAlpha.alphaR = 0.05;
    else
        if testAlpha.alphaR < 0 || testAlpha.alphaR > 1
            disp('--detectSubResFeatures2D_V2: testAlpha.alphaR should be between 0 and 1!');
            errFlag = 1;
        end
    end
    if ~isfield(testAlpha,'alphaA')
        testAlpha.alphaA = 0.05;
    else
        if testAlpha.alphaA < 0 || testAlpha.alphaA > 1
            disp('--detectSubResFeatures2D_V2: testAlpha.alphaA should be between 0 and 1!');
            errFlag = 1;
        end
    end
    if ~isfield(testAlpha,'alphaD')
        testAlpha.alphaD = 0.05;
    else
        if testAlpha.alphaD < 0 || testAlpha.alphaD > 1
            disp('--detectSubResFeatures2D_V2: testAlpha.alphaD should be between 0 and 1!');
            errFlag = 1;
        end
    end
    if ~isfield(testAlpha,'alphaF')
        testAlpha.alphaF = 0.05;
    else
        if testAlpha.alphaF < 0 || testAlpha.alphaF > 1
            disp('--detectSubResFeatures2D_V2: testAlpha.alphaF should be between 0 and 1!');
            errFlag = 1;
        end
    end
    
end

%check visualization option
if nargin < 5 || isempty(visual)
    visual = 0;
else
    if visual ~= 0 && visual ~= 1
        disp('--detectSubResFeatures2D_V2: Variable "visual" should be 0 or 1!');
        errFlag = 1;
    end
end

%check whether to do MMF
if nargin < 6 || isempty(doMMF)
    doMMF = 1;
else
    if doMMF ~= 0 && doMMF ~= 1
        disp('--detectSubResFeatures2D_V2: Variable "doMMF" should be 0 or 1!');
        errFlag = 1;
    end
end

%check the bit depth
if nargin < 7 || isempty(bitDepth)
    bitDepth = 14;
else
    if bitDepth <= 0 || bitDepth-floor(bitDepth) ~= 0
        disp('--detectSubResFeatures2D_V2: Variable "bitDepth" should be a positive integer!');
    end
end

%check whether results are to be saved
if nargin < 8 || isempty(saveResults)
    saveResults = 0;
else
    if saveResults ~= 0 && saveResults ~= 1
        disp('--detectSubResFeatures2D_V2: Variable "saveResults" should be 0 or 1!');
    end
end

%check whether back ground noise sigma is input or whether it should be
%calculated on the fly
if nargin < 9 || isempty(bgNoiseSigma)
    bgNoiseSigma = 0;
    estimateBgNoise = 1;
else
    estimateBgNoise = 0;
end

%exit if there are problems with input data
if errFlag
    disp('--detectSubResFeatures2D_V2: Please fix input data!');
    return
end

%get number of pixels in each direction (in image coordinate system)
% [numPixelsY,numPixelsX] = size(image);
% [numPixelsY,numPixelsX,numPixelsZ] = size(image);
[numPixelsX,numPixelsY,numPixelsZ] = size(image);

%extract test alpha values from input
alphaR = testAlpha.alphaR;
alphaA = testAlpha.alphaA;
alphaD = testAlpha.alphaD;
alphaF = testAlpha.alphaF;

%Divide image by bit depth, to normalize it between 0 and 1
image = double(image)/(2^bitDepth-1);

%get background intensity information from cands
bgAmp = vertcat(cands.IBkg);
status = vertcat(cands.status);
bgAmp = bgAmp(status==1);
bgAmpAve = mean(bgAmp);

%% Determine signal overlap

%determine which signals are overlapping, in which case they must
%be fitted together later on
% [clusters,errFlag] = findOverlapPSFs2D(cands,numPixelsX,numPixelsY,psfSigma);
% if errFlag
%     disp('--detectSubResFeatures2D_V2: Could not place signals in clusters!');
%     return
% end
%[clusters,errFlag] = findOverlapPSFs3D(cands,numPixelsX,numPixelsY,numPixelsZ,psfSigma);
[clusters,errFlag] = findOverlapPSFs3D_new_reCluster(cands,numPixelsX,numPixelsY,numPixelsZ,psfSigma);
if errFlag
    disp('--detectSubResFeatures3D: Could not place signals in clusters!');
    return
end

%% Mixture Model Fitting

%set optimization options
options = optimset('Jacobian','on',...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

alphaR_def = alphaR;

% recluster while loop
reCluster = 1;
while reCluster
    % set initially to no reCluster
    reCluster = 0;
    % results matrix (x,y,z,A,dx,dy,dz,dA)
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
        bgAmpT = bgAmpAve;
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
            [solutionT,residualsT,jacMatT] = fitNGaussians3D_mexCode_fitFun(options,x0, lb, ub, imageC, clusterPixels, psfSigma, [numPixelsX numPixelsY numPixelsZ], clusterPixels1DIdx);
            
            jacMatT = full(jacMatT);
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
                    
                    disp('spot added')
                    
                    clusters(iCluster).maximaPos = [solutionT(1:4:end-1) solutionT(2:4:end-1) solutionT(3:4:end-1)];
                    clusters(iCluster).maximaAmp = solutionT(4:4:end-1);
                    clusters(iCluster).numMaxima = numMaximaT;
                    % re-cluster
                    clusters_tmp = findOverlapPSFs3D_new_reCluster(clusters,numPixelsX,numPixelsY,numPixelsZ,psfSigma);
                    % check number of clusters
                    if length(clusters_tmp) < length(clusters)% if so then keep the new cluster struct and re-start the fitting
                        clusters = clusters_tmp;
                        reCluster = 1;
                        fit = 0;
                    else % else just fit current cluster and get new pixels to fit (if any), and re-fit
                        clusters_tmp = findOverlapPSFs3D_new_reCluster(clusters(iCluster),numPixelsX,numPixelsY,numPixelsZ,psfSigma);
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
                            [solutionT,residualsT,jacMatT] = fitNGaussians3D_mexCode_fitFun(options,solutionT, lb, ub, imageC, clusterPixels, psfSigma, [numPixelsX numPixelsY numPixelsZ], clusterPixels1DIdx);
                            jacMatT = full(jacMatT);
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
                jacMat = jacMatT;
                
                %calculate the parameters' variance-covariance matrix and get their
                %uncertainties
                varCovMat = residVar * inv(jacMat'*jacMat);
                standDevVec = sqrt(diag(varCovMat));
                
                %if nothing weird happened in the fit...
                if all(isreal(standDevVec))
                    
                    %extract estimate and std of background intensity and
                    %remove from vectors
                    bgAmp = [solution(end) standDevVec(end)];
                    solution = solution(1:end-1);
                    standDevVec = standDevVec(1:end-1);
                    
                    %                 %reshape 3nx1 vectors "solution" and "standDevVec" into nx3 matrices
                    %                 solution = reshape(solution,3,numMaxima);
                    %                 solution = solution';
                    %                 standDevVec = reshape(standDevVec,3,numMaxima);
                    %                 standDevVec = standDevVec';
                    
                    %reshape 4nx1 vectors "solution" and "standDevVec" into nx3 matrices
                    solution = reshape(solution,4,numMaxima);
                    solution = solution';
                    standDevVec = reshape(standDevVec,4,numMaxima);
                    standDevVec = standDevVec';
                    
                    %                 %extract feature positions and amplitudes and their uncertainties
                    %                 maximaPos = [solution(:,1:2) standDevVec(:,1:2)];
                    %                 maximaAmp = [solution(:,3) standDevVec(:,3)];
                    
                    %extract feature positions and amplitudes and their uncertainties
                    maximaPos = [solution(:,1:3) standDevVec(:,1:3)];
                    maximaAmp = [solution(:,4) standDevVec(:,4)];
                    
                    %if attempting iterative mixture-model fitting
                    if doMMF
                        
                        %add a kernel positioned at the pixel with maximum residual
                        %for fit with numMaxima + 1 Gaussians
                        numMaximaT = numMaxima + 1; %update number of maxima
                        maximaAmpT = [maximaAmp(:,1); mean(maximaAmp(:,1))]; %signal amplitude
                        indxMaxRes = find(residuals==max(residuals)); %index of pixel with maximum residuals
                        coord = clusterPixels(indxMaxRes(1),:); %position of new kernel (using indxMaxRes(1) avoids crashing if it so happens that two pixels have the exact same value which happens to be the maximum value - it has happened!!!)
                        %                     maximaPosT = [maximaPos(:,1:2); coord]; %signal positions
                        maximaPosT = [maximaPos(:,1:3); coord]; %signal positions
                        bgAmpT = bgAmp(1); %background amplitude
                        
                    else %if no iterative mixture model fitting
                        
                        fit = 0;
                        
                    end
                    
                else %if something weird did happen in the fit
                    
                    fit = 0; %stop fitting this cluster
                    keepCluster = 0; %mark it to be discarded
                    
                end %(if all(isreal(standDevVec)))
                
            end %(if fit)
            
        end %(while fit)
        
        % if reClustering break here
        if reCluster
            break
        end
        
        %if nothing weird happened in the fit and this cluster is to be
        %retained
        if keepCluster
            
            %% check amplitudes and make sure that they are all significant
            
            %1-sided t-test: H0: T=0, H1: T>0
            %calculate test statistic (t-distributed)
            testStat = maximaAmp(:,1)./sqrt((maximaAmp(:,2).^2+residVar));
            
            %get p-value
            pValue = 1-tcdf(testStat,numDegFree);
            
            %find largest p-value and decide whether to remove one kernel,
            %repeat the fit and test again for amplitude
            pValueMax = max(pValue);
            testAmp = pValueMax > alphaA;
            
            %if any of the amplitudes is not significant
            while testAmp
                
                %remove a kernel and refit only if there is originally more
                %than one kernel
                %otherwise, discard the whole cluster
                if numMaxima > 1
                    
                    %determine which kernels to keep
                    indxBad = find(pValue == pValueMax);
                    indxBad = indxBad(1);
                    indx = setdiff(1:numMaxima,indxBad);
                    
                    %remove the information of the kernel to be discarded
                    maximaPos = maximaPos(indx,:);
                    maximaAmp = maximaAmp(indx,:);
                    numMaxima = numMaxima - 1;
                    
                    %collect initial guesses and lower and upper bounds for fit
                    %                 [x0,lb,ub] = mmfInitGuessLowerUpperBounds(maximaPos(:,1:2),...
                    %                     maximaAmp(:,1),bgAmp(1),psfSigma,clusterPixels,1);
                    
                    %clusterPixels = clusterPixels_new;
                    [x0,lb,ub] = mmfInitGuessLowerUpperBounds3D(maximaPos(:,1:3),...
                        maximaAmp(:,1),bgAmp(1),psfSigma,clusterPixels,1);
                    
                    %calculate number of degrees of freedom in system
                    %                 numDegFree = size(clusterPixels,1) - 3*numMaxima - 1;
                    
                    numDegFree = size(clusterPixels,1) - 4*numMaxima - 1;
                    
                    %determine feature positions and amplitudes and estimate background
                    %intensity, using nonlinear least squares data fitting
                    %also get residuals and Jacobian matrix to calculate variance-covariance
                    %of estimated parameters
                    %                 [solution,~,residuals,~,~,~,jacMat] = ...
                    %                     lsqnonlin(@fitNGaussians2D,x0,lb,ub,options,imageC,...
                    %                     clusterPixels,psfSigma);
                    
                    %[solution,~,residuals,~,~,~,jacMat] = ...
                    %    lsqnonlin(@fitNGaussians3D_mexCode_mex,x0,lb,ub,options,imageC,...
                    %    clusterPixels,psfSigma);
                    
                    [solution,residuals,jacMat] = fitNGaussians3D_mexCode_fitFun(options, x0, lb, ub, imageC, clusterPixels, psfSigma, [numPixelsX,numPixelsY,numPixelsZ], clusterPixels1DIdx);
                    jacMat = full(jacMat);
                    residuals = -residuals; %minus sign so that residuals = real image - model image
                    residVar = sum(residuals.^2)/numDegFree;
                    
                    %calculate the parameters' variance-covariance matrix and get their
                    %uncertainties
                    varCovMat = residVar * inv(jacMat'*jacMat);
                    standDevVec = sqrt(diag(varCovMat));
                    
                    %if nothing weird happened in the fit...
                    if all(isreal(standDevVec))
                        
                        % test for merging clusters
                        clusters(iCluster).maximaPos = [solution(1:4:end-1) solution(2:4:end-1) solution(3:4:end-1)];
                        clusters(iCluster).maximaAmp = solution(4:4:end-1);
                        clusters(iCluster).numMaxima = numMaxima;
                        % re-cluster
                        clusters_tmp = findOverlapPSFs3D_new_reCluster(clusters,numPixelsX,numPixelsY,numPixelsZ,psfSigma);
                        
                        if length(clusters_tmp) < length(clusters)
                            clusters = clusters_tmp;
                            reCluster = 1;
                            break
                        else % else get new list of pixels
                            clusters_tmp = findOverlapPSFs3D_new_reCluster(clusters(iCluster),numPixelsX,numPixelsY,numPixelsZ,psfSigma);
                            if length(clusters_tmp) > 1
                                tmp = clusters_tmp;
                                clear clusters_tmp
                                tmp1 = catStruct(1,'tmp.pixels');
                                clusters_tmp.pixels = tmp1;
                                clear tmp1 tmp
                            end
                            clusterPixels = clusters_tmp.pixels(:,1:3);
                            clusterPixels1DIdx = clusters_tmp.pixels(:,4);
                            imageC = image(clusterPixels1DIdx);
                        end
                        
                        %extract estimate and std of background intensity and
                        %remove from vectors
                        bgAmp = [solution(end) standDevVec(end)];
                        solution = solution(1:end-1);
                        standDevVec = standDevVec(1:end-1);
                        
                        %                     %reshape 3nx1 vectors "solution" and "standDevVec" into nx3 matrices
                        %                     solution = reshape(solution,3,numMaxima);
                        %                     solution = solution';
                        %                     standDevVec = reshape(standDevVec,3,numMaxima);
                        %                     standDevVec = standDevVec';
                        
                        %reshape 4nx1 vectors "solution" and "standDevVec" into nx3 matrices
                        solution = reshape(solution,4,numMaxima);
                        solution = solution';
                        standDevVec = reshape(standDevVec,4,numMaxima);
                        standDevVec = standDevVec';
                        
                        %                     %extract feature positions and amplitudes and their uncertainties
                        %                     maximaPos = [solution(:,1:2) standDevVec(:,1:2)];
                        %                     maximaAmp = [solution(:,3) standDevVec(:,3)];
                        
                        %extract feature positions and amplitudes and their uncertainties
                        maximaPos = [solution(:,1:3) standDevVec(:,1:3)];
                        maximaAmp = [solution(:,4) standDevVec(:,4)];
                        
                        %check amplitudes and make sure that they are all significant
                        
                        %1-sided t-test: H0: T=0, H1: T>0
                        %calculate test statistic (t-distributed)
                        testStat = maximaAmp(:,1)./sqrt((maximaAmp(:,2).^2+residVar));
                        
                        %get p-value
                        pValue = 1-tcdf(testStat,numDegFree);
                        
                        %find largest p-value and decide whether to remove one kernel,
                        %repeat the fit and test again for amplitude
                        pValueMax = max(pValue);
                        testAmp = pValueMax > alphaA;
                        
                    else %(if all(isreal(standDevVec)))
                        
                        testAmp = 0; %stop fitting this cluster
                        keepCluster = 0; %mark it to be discarded
                        
                    end %(if all(isreal(standDevVec)))
                    
                else %(if numMaxima > 1)
                    
                    testAmp = 0; %stop fitting this cluster
                    keepCluster = 0; %mark it to be discarded
                    
                end %(if numMaxima > 1)
                
            end %(while testAmp)
            
            % if reClustering break here
            if reCluster
                break
            end
            
            %% if there is more than one maximum, test for sigificance of distances between maxima
            
            if numMaxima > 1
                
                %get p-values of distances between maxima
                %             pValue = mmfDistPV(maximaPos,varCovMat,numMaxima,numDegFree);
                pValue = mmfDistPV3D(maximaPos,varCovMat,numMaxima,numDegFree);
                
                %find largest p-value and decide whether to remove one kernel,
                %repeat the fit and test again for distances
                pValueMax = max(pValue(:));
                testDist = pValueMax > alphaD;
                
                %if any of the distances is not significant
                while testDist
                    
                    %find pair with maximum p-value
                    [indx1,indx2] = find(pValue==pValueMax);
                    indx1 = indx1(1);
                    indx2 = indx2(1);
                    
                    %out of this pair, identify maximum with smaller amplitude
                    ampBoth = maximaAmp([indx1 indx2],1);
                    if ampBoth(1) < ampBoth(2)
                        indxBad = indx1;
                    else
                        indxBad = indx2;
                    end
                    
                    %determine which kernels to keep
                    indx = setdiff(1:numMaxima,indxBad);
                    
                    %remove the information of the kernel to be discarded
                    maximaPos = maximaPos(indx,:);
                    maximaAmp = maximaAmp(indx,:);
                    numMaxima = numMaxima - 1;
                    
                    %collect initial guesses and lower and upper bounds for fit
                    %                 [x0,lb,ub] = mmfInitGuessLowerUpperBounds(maximaPos(:,1:2),...
                    %                     maximaAmp(:,1),bgAmp(1),psfSigma,clusterPixels,1);
                    
                    %clusterPixels = clusterPixels_new;
                    [x0,lb,ub] = mmfInitGuessLowerUpperBounds3D(maximaPos(:,1:3),...
                        maximaAmp(:,1),bgAmp(1),psfSigma,clusterPixels,1);
                    
                    
                    %calculate number of degrees of freedom in system
                    %                 numDegFree = size(clusterPixels,1) - 3*numMaxima - 1;
                    
                    numDegFree = size(clusterPixels,1) - 4*numMaxima - 1;
                    
                    %determine feature positions and amplitudes and estimate background
                    %intensity, using nonlinear least squares data fitting
                    %also get residuals and Jacobian matrix to calculate variance-covariance
                    %of estimated parameters
                    [solution,residuals,jacMat] = fitNGaussians3D_mexCode_fitFun(options, x0, lb, ub, imageC, clusterPixels, psfSigma, [numPixelsX,numPixelsY,numPixelsZ], clusterPixels1DIdx);
                    jacMat = full(jacMat);
                    residuals = -residuals; %minus sign so that residuals = real image - model image
                    residVar = sum(residuals.^2)/numDegFree;
                    
                    %calculate the parameters' variance-covariance matrix and get their
                    %uncertainties
                    varCovMat = residVar * inv(jacMat'*jacMat);
                    standDevVec = sqrt(diag(varCovMat));
                    
                    %if nothing weird happened in the fit...
                    if all(isreal(standDevVec))
                        
                        % test for merging clusters
                        clusters(iCluster).maximaPos = [solution(1:4:end-1) solution(2:4:end-1) solution(3:4:end-1)];
                        clusters(iCluster).maximaAmp = solution(4:4:end-1);
                        clusters(iCluster).numMaxima = numMaxima;
                        % re-cluster
                        clusters_tmp = findOverlapPSFs3D_new_reCluster(clusters,numPixelsX,numPixelsY,numPixelsZ,psfSigma);
                        
                        if length(clusters_tmp) < length(clusters)
                            clusters = clusters_tmp;
                            reCluster = 1;
                            break
                        else % else get new list of pixels
                            clusters_tmp = findOverlapPSFs3D_new_reCluster(clusters(iCluster),numPixelsX,numPixelsY,numPixelsZ,psfSigma);
                            if length(clusters_tmp) > 1
                                tmp = clusters_tmp;
                                clear clusters_tmp
                                tmp1 = catStruct(1,'tmp.pixels');
                                clusters_tmp.pixels = tmp1;
                                clear tmp1 tmp
                            end
                            clusterPixels = clusters_tmp.pixels(:,1:3);
                            clusterPixels1DIdx = clusters_tmp.pixels(:,4);
                            imageC = image(clusterPixels1DIdx);
                        end
                        
                        
                        
                        %extract estimate and std of background intensity and
                        %remove from vectors
                        bgAmp = [solution(end) standDevVec(end)];
                        solution = solution(1:end-1);
                        standDevVec = standDevVec(1:end-1);
                        
                        %                     %reshape 3nx1 vectors "solution" and "standDevVec" into nx3 matrices
                        %                     solution = reshape(solution,3,numMaxima);
                        %                     solution = solution';
                        %                     standDevVec = reshape(standDevVec,3,numMaxima);
                        %                     standDevVec = standDevVec';
                        
                        %reshape 3nx1 vectors "solution" and "standDevVec" into nx3 matrices
                        solution = reshape(solution,4,numMaxima);
                        solution = solution';
                        standDevVec = reshape(standDevVec,4,numMaxima);
                        standDevVec = standDevVec';
                        
                        %extract feature positions and amplitudes and their uncertainties
                        %                     maximaPos = [solution(:,1:2) standDevVec(:,1:2)];
                        %                     maximaAmp = [solution(:,3) standDevVec(:,3)];
                        
                        maximaPos = [solution(:,1:3) standDevVec(:,1:3)];
                        maximaAmp = [solution(:,4) standDevVec(:,4)];
                        
                        
                        %if there is more than one maximum, test for
                        %sigificance of distances between maxima
                        if numMaxima > 1
                            
                            %get p-values of distances between maxima
                            %                         pValue = mmfDistPV(maximaPos,varCovMat,numMaxima,numDegFree);
                            pValue = mmfDistPV3D(maximaPos,varCovMat,numMaxima,numDegFree);
                            
                            %find largest p-value and decide whether to remove one kernel,
                            %repeat the fit and test again for distances
                            pValueMax = max(pValue(:));
                            testDist = pValueMax > alphaD;
                            
                        else
                            
                            testDist = 0;
                            
                        end %(if numMaxima > 1)
                        
                    else %(if all(isreal(standDevVec)))
                        
                        testDist = 0; %stop fitting this cluster
                        keepCluster = 0; %mark it to be discarded
                        
                    end %(if all(isreal(standDevVec)))
                    
                end %(while testDist)
                
                % if reClustering break here
                if reCluster
                    break
                end
                
            end %(if numMaxima > 1)
            
        else %(if all(isreal(standDevVec)))
            
            keepCluster = 0; %mark this cluster to be discarded
            
        end %(if all(isreal(standDevVec)))
        
        
        %% estimate of the parameters and the residuals from the fit
        if keepCluster
            
            %         %collect initial guesses and lower and upper bounds for fit
            %         [x0,lb,ub] = mmfInitGuessLowerUpperBounds(maximaPos(:,1:2),...
            %             maximaAmp(:,1),bgAmp(1),psfSigma,clusterPixels,1);
            %collect initial guesses and lower and upper bounds for fit
            %clusterPixels = clusterPixels_new;
            [x0,lb,ub] = mmfInitGuessLowerUpperBounds3D(maximaPos(:,1:3),...
                maximaAmp(:,1),bgAmp(1),psfSigma,clusterPixels,1);
            
            
            %calculate number of degrees of freedom in system
            %         numDegFree = size(clusterPixels,1)-3*numMaxima-1;
            
            numDegFree = size(clusterPixels,1)-4*numMaxima-1;
            
            %determine feature positions and amplitudes and estimate background
            %intensity, using nonlinear least squares data fitting
            %also get residuals and Jacobian matrix to calculate variance-covariance
            %of estimated parameters
            
            [solution,residuals,jacMat] = fitNGaussians3D_mexCode_fitFun(options, x0, lb, ub, imageC, clusterPixels, psfSigma, [numPixelsX,numPixelsY,numPixelsZ], clusterPixels1DIdx);
            
            % test for merging clusters
            clusters(iCluster).maximaPos = [solution(1:4:end-1) solution(2:4:end-1) solution(3:4:end-1)];
            clusters(iCluster).maximaAmp = solution(4:4:end-1);
            clusters(iCluster).numMaxima = numMaxima;
            % re-cluster
            clusters_tmp = findOverlapPSFs3D_new_reCluster(clusters,numPixelsX,numPixelsY,numPixelsZ,psfSigma);
            
            if length(clusters_tmp) < length(clusters)
                clusters = clusters_tmp;
                reCluster = 1;
                break
            end
            
            jacMat = full(jacMat);
            residuals = -residuals; %minus sign so that residuals = real image - model image
            residVar = sum(residuals.^2)/numDegFree;
            
            %calculate the parameters' variance-covariance matrix and get their
            %uncertainties
            varCovMat = residVar * inv(jacMat'*jacMat);
            standDevVec = sqrt(diag(varCovMat));
            
            %extract estimate and std of background intensity and
            %remove from vectors
            bgAmp = [solution(end) standDevVec(end)];
            solution = solution(1:end-1);
            standDevVec = standDevVec(1:end-1);
            
            %         %reshape 3nx1 vectors "solution" and "standDevVec" into nx3 matrices
            %         solution = reshape(solution,3,numMaxima);
            %         solution = solution';
            %         standDevVec = reshape(standDevVec,3,numMaxima);
            %         standDevVec = standDevVec';
            
            %reshape 4nx1 vectors "solution" and "standDevVec" into nx4 matrices
            solution = reshape(solution,4,numMaxima);
            solution = solution';
            standDevVec = reshape(standDevVec,4,numMaxima);
            standDevVec = standDevVec';
            
            %         %extract feature positions and amplitudes and their uncertainties
            %         maximaPos = [solution(:,1:2) standDevVec(:,1:2)];
            %         maximaAmp = [solution(:,3) standDevVec(:,3)];
            
            %extract feature positions and amplitudes and their uncertainties
            maximaPos = [solution(:,1:3) standDevVec(:,1:3)];
            maximaAmp = [solution(:,4) standDevVec(:,4)];
            
            %             %store solution in clustersMMF
            %             clustersMMF(i).position = maximaPos;
            %             clustersMMF(i).amplitude = maximaAmp;
            %             clustersMMF(i).bgAmp = bgAmp;
            %             clustersMMF(i).numDegFree = numDegFree;
            %             clustersMMF(i).residuals = residuals;
            
            %estimate background noise variance if not supplied
            bgNoiseVar = bgNoiseSigma^2;
            %go over cluster and check that its residuals are
            %comparable to the background noise
            testStat = (sum(residuals.^2)/numDegFree)/bgNoiseVar;
            %get p-value of test statistic
            pValue = 1 - fcdf(testStat,numDegFree,length(imageC)-1);
            %compare p-value to alpha
            %1-sided F-test: H0: F=1, H1: F>1
            if pValue < alphaF %if p-value is smaller than alpha, reject overall fit
                keepCluster = 0;
            end
            
            
        end %(if keepCluster)
        
        % discard cluster if not keeping it and reset cluster count
        if ~keepCluster
            clusters(iCluster) = [];
            iCluster = iCluster - 1;
            continue
        end
        
        % add results to spots matrix
        
        spots = [spots; maximaPos(:,1:3) maximaAmp(:,1) maximaPos(:,4:6) maximaAmp(:,2)];
        
    end %while iCluster < length(clusters)
end % while reClusrer

%% results

if ~isempty(spots)
    
    %store information in structure "detectedFeatures"
    detectedFeatures.xCoord = [spots(:,1) spots(:,5)];
    detectedFeatures.yCoord = [spots(:,2) spots(:,6)];
    detectedFeatures.zCoord = [spots(:,3) spots(:,7)];
    detectedFeatures.amp = [spots(:,4) spots(:,8)];
    
    
else
    
    %store empty structures
    detectedFeatures = struct('xCoord',zeros(0,2),'yCoord',zeros(0,2),'zCoord',zeros(0,2),'amp',zeros(0,2));
    
end

%save output if requested
if saveResults
    save('detectedFeatures','detectedFeatures','clustersMMF');
end


%% Visualization

if visual
    
    %     %make 3 layers out of original image (normalized)
    %     imageNorm = image/max(image(:));
    %     imageN3 = repmat(imageNorm,[1 1 3]);
    %
    %     %place zeros in pixels of maxima from cands
    %     for i=1:length(clusters)
    %         for j=1:3
    %             imageN3(clusters(i).maximaPos(:,3)+(j-1)*numPixelsX*numPixelsY)=0;
    %         end
    %     end
    %
    %     %place zeros in pixels of maxima from mixture-model fitting
    %     for i=1:length(clustersMMF)
    %         pos = (round(clustersMMF(i).position(:,1)-1))*numPixelsY ...
    %             + round(clustersMMF(i).position(:,2));
    %         for j=1:3
    %             imageN3(pos+(j-1)*numPixelsX*numPixelsY)=0;
    %         end
    %     end
    %
    %     %label maxima from cands in blue
    %     for i=1:length(clusters)
    %         imageN3(clusters(i).maximaPos(:,3)+2*numPixelsX*numPixelsY)=1;
    %     end
    %
    %     %label maxima from mixture-model fitting in red
    %     %a maximum from mixture-model fitting that falls in the same pixel
    %     %as that from cands will appear in magenta
    %     for i=1:length(clustersMMF)
    %         pos = (round(clustersMMF(i).position(:,1)-1))*numPixelsY ...
    %             + round(clustersMMF(i).position(:,2));
    %         imageN3(pos)=1;
    %     end
    %
    %     %plot image
    %     imtool(imageN3);
    
end

%%%%% ~~ the end ~~ %%%%%

