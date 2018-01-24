function dataStruct = makiFastMMF3D(dataStruct,varargin)
% finds spots in 3D using fast MMF spotFit.m, based on detectSpots3D.m,
% uses splitModes.m for localMaxima
% also extimes psf of movie in 3D

%% Input + initialization

%check whether correct number of input arguments was used
if nargin < 1
    disp('--detectSpots: Incorrect number of input arguments!');
    return
end

% extra param
param = struct('numSigmaIter',[],'integWindow',[]);
fn = fieldnames(param);
%I loop through the variables names given in Parameters.
for i = 1:2:length(varargin)%(variables names are odd value of varargin
    if validatestring(varargin{i},fn)
        param.(varargin{i}) = varargin{i+1};
    end%If i find the value i replace the default value.
end

%get movie parameters
imageName = fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName);

% get number of frames
nFrames = dataStruct.dataProperties.movieSize(4);

% use deconvolved image or not
if isfield(dataStruct.dataProperties,'decon') && ~isempty(dataStruct.dataProperties.decon)
    decon = dataStruct.dataProperties.decon;
else
    decon = 0;
end

% get cropping if any
if isfield(dataStruct.dataProperties,'crop')
    crop = dataStruct.dataProperties.crop;
else
    crop = [];
end

%get initial guess of PSF sigma
filterPrm = dataStruct.dataProperties.FILTERPRM;

%get camera bit depth
if ~isfield(dataStruct.dataProperties,'bitDepth') || isempty(dataStruct.dataProperties.bitDepth)
    bitDepth = 16;
else
    bitDepth = dataStruct.dataProperties.bitDepth;
end

%check whether to estimate PSF sigma from the data, take first from
%varagrin, and then from dataProperties, and then a default value
if isempty(param.numSigmaIter)
    if ~isfield(dataStruct.dataProperties,'numSigmaIter') || isempty(dataStruct.dataProperties.numSigmaIter)
        numSigmaIter = 10;
    else
        numSigmaIter = dataStruct.dataProperties.numSigmaIter;
    end
else
    numSigmaIter = param.numSigmaIter;
end

%get integration time window, take first from
%varagrin, and then from dataProperties, and then a default value
if isempty(param.integWindow)
    if ~isfield(dataStruct.dataProperties,'integWindow')
        integWindow = 0;
    else
        integWindow = dataStruct.dataProperties.integWindow;
    end
else
    integWindow = param.integWindow;
end
numIntegWindow = length(integWindow);

%initialize some variables
emptyFrames = [];
framesFailedLocMax = [];
framesFailedMMF = [];

%turn warnings off
warningState = warning('off','all');

%% General image information

numImagesInteg = repmat(nFrames,1,numIntegWindow) - 2 * integWindow; %integrated images

% read image from .mat file
try
    movie = readOMEMatFile(imageName,1:nFrames,1,decon,crop);
    [imageSizeX,imageSizeY,imageSizeZ,~] = size(movie);
catch %#ok
    error('Cannot read .mat image file');
end

% divide by bit depth
movie = movie ./ (2^bitDepth - 1);

%% Local maxima detection

%initialize output structure
localMaxima = repmat(struct('cands',[]),nFrames,1);

% crate signal filter
signalFilter = GaussMask3D(filterPrm(1:3),filterPrm(4:6),[],1,[],[],1);

% create background filter
backgroundFilterParms = filterPrm(1:3) * 14;
backgroundFilterParms(4:6) = roundOddOrEven(backgroundFilterParms,'odd','inf');
backgroundFilter = GaussMask3D(backgroundFilterParms(1:3),backgroundFilterParms(4:6),[],1,[],[],1);

% make separated noise mask
% noiseMask = {ones(filterPrm(4),1,1)./filterPrm(4),ones(1,filterPrm(5),1)./filterPrm(5),ones(1,1,filterPrm(6))./filterPrm(6)};

% save index of each frame for image integration
imageIndx = 1:nFrames;

for iWindow = 1 : numIntegWindow
    
    %initialize progress text
    progressText(0,['Detecting local maxima with integration window = ' num2str(integWindow(iWindow))]);
    
    for iImage = 1 : numImagesInteg(iWindow)
        
        % get frame
        imageRaw = movie(:,:,:,imageIndx((1:1+2*integWindow(iWindow))+iImage-1));
        
        %replace zeros with NaNs
        %zeros result from cropping that leads to curved boundaries
        imageRaw(imageRaw==0) = NaN;
        
        %integrate images
        imageInteg = nanmean(imageRaw,4);
        
        %filter integrated image
        imageIntegF = fastGauss3D(imageInteg,[],filterPrm(4:6),1,signalFilter,any(isnan(imageInteg(:))));
        
        % get first guess at background
        background = fastGauss3D(imageIntegF,[],backgroundFilterParms(4:6),1,backgroundFilter,any(isnan(imageIntegF(:))));
        
        % mask signal pixels, then recalculate background, define a signal
        % pixel as a pixel in the raw image that is greater than 4* the
        % background
        sigMask = imageInteg > 4*background;
        
        % remove some of the spurious hits. Use 2d mask for speed and
        % memory (3d strel won't work on binary image with imopen)
        sigMask = imopen(sigMask,strel('disk',1));
        
        % fill in binary mask holes here
        sigMask = imfill(sigMask,'holes');
        
        % create raw mask
        filtMsk = imageIntegF.*~sigMask;
        clear sigMask
        
        % recalculate background
        filtMsk(filtMsk==0) = NaN;
        background = fastGauss3D(filtMsk,[],backgroundFilterParms(4:6),1,backgroundFilter,any(isnan(filtMsk(:))));
        clear filtMsk
        
        % put zeros back
        background(isnan(background)) = 0;
        imageIntegF(isnan(imageIntegF)) = 0;
        % imageInteg(isnan(imageInteg)) = 0;
        
        % amplitude is filtered - background
        amplitude = imageIntegF - background;
        clear background
        
        % noise is local average (averaged over filter support) of squared
        % residuals of raw-filtered image
        % noise = (imageInteg - imageIntegF).^2;
        % noise = fastGauss3D(noise,[],filterPrm(4:6),1,noiseMask);
        clear imageInteg imageIntegF
        
        try
            
            % get local maxima from the amplitude
            locMax = loc_max3Df(amplitude);
            localMax1DIndx = sub2ind([imageSizeX imageSizeY imageSizeZ],locMax(:,1),locMax(:,2),locMax(:,3));
            
            % get signal to noise of local maxima, amplitude/sqrt(nse/amp) - poisson
            % s2n = amplitude(localMax1DIndx) ./ sqrt(noise(localMax1DIndx) ./ amplitude(localMax1DIndx));
            amp = amplitude(localMax1DIndx);
            clear amplitude localMax1DIndx
            
            % set s2n with complex values to zero
            % realVals = arrayfun(@(x) isreal(x(1)),s2n);
            % s2n(~realVals) = 0;
            % clear realVals
            
            % get amp cutOff using unimodal thresholding, only use value
            % above zero
            cutVal = splitModes(amp(amp>0));
            
            % get indexes of spots above cutOff
            passIdx = amp > cutVal;
            
            % keep these spots
            cands = locMax(passIdx,:);
            % clear amp passIdx locMax
            clear amp locMax
            
            %add the cands of the current image to the rest - this is done
            %for the raw images, not the integrated ones
            localMaxima(iImage+integWindow(iWindow)).cands = ...
                [localMaxima(iImage+integWindow(iWindow)).cands; cands];
            
        catch %#ok<CTCH> if locMax fails
            framesFailedLocMax = [framesFailedLocMax; iImage]; %#ok<AGROW>
        end
        
        %display progress
        progressText(iImage/numImagesInteg(iWindow),['Detecting local maxima with integration window = ' num2str(integWindow(iWindow))]);
        
    end %(for iImage = 1 : numImagesInteg(iWindow))
    
    %assign local maxima for frames left out due to time integration
    for iImage = 1 : integWindow(iWindow)
        localMaxima(iImage).cands = [localMaxima(iImage).cands; ...
            localMaxima(integWindow(iWindow)+1).cands];
    end
    for iImage = nFrames-integWindow(iWindow)+1 : nFrames
        localMaxima(iImage).cands = [localMaxima(iImage).cands; ...
            localMaxima(end-integWindow(iWindow)).cands];
    end
    
end %(for iWindow = 1 : numIntegWindow)

% calculate s2n cutoff, don't include zero value as they give artificially
% low cutoffs
%s2n = catStruct(1,'localMaxima.cands(:,4)');
%cutVal = splitModes(s2n(s2n>0));
%clear s2n

% get only unique failed frames if any
framesFailedLocMax = unique(framesFailedLocMax);

%go over all frames, perform thresholing, remove redundant cands, and register empty frames
progressText(0,'Removing redundant local maxima');
for iImage = 1 : nFrames
    
    %get the cands of this frame
    candsCurrent = localMaxima(iImage).cands;
    
    %if there are no cands, register that this is an empty frame
    if isempty(candsCurrent)
        
        emptyFrames = [emptyFrames; iImage]; %#ok<AGROW>
        
    else
        
        % perform thresholding
        %candsCurrent = candsCurrent(candsCurrent(:,4) > cutVal,1:3);
        
        %find the unique local maxima positions
        candsCurrent = unique(candsCurrent,'rows');
        
        %keep only these unique cands
        localMaxima(iImage).cands = candsCurrent;
        
    end
    
    %display progress
    progressText(iImage/nFrames,'Removing redundant local maxima');
    
end

%make a list of images that have local maxima
goodImages = setxor(1:nFrames,emptyFrames)';

% get psf sigma from filterPrm
psfSigma = [filterPrm(1) filterPrm(3)];

%% PSF sigma estimation

if numSigmaIter
    
    %specify which parameters to fit for sigma estimation
    %         fitParameters = [{'X1'} {'X2'} {'A'} {'Sxy'} {'B'}];
    fitParameters = [{'X1'} {'X2'} {'X3'} {'A'} {'Sxy'} {'S3'} {'B'}]; % 3d approx
    
    %store original input sigma
    psfSigmaIn = psfSigma;
    
    %give a dummy value for psfSigma0 and acceptCalc to start while loop
    %psfSigma0 = 0;
    psfSigma0 = [0 0];
    acceptCalc = 1;
    
    %initialize variable counting number of iterations
    numIter = 0;
    
    %iterate as long as estimated sigma is larger than initial sigma
    %while numIter <= numSigmaIter && acceptCalc && ((psfSigma-psfSigma0)/psfSigma0 > 0.05)
    while numIter <= numSigmaIter && acceptCalc && any((abs(psfSigma - psfSigma0) ./ psfSigma0) > [0.05 0.05]) % fit both sigma-xy and sigma-z
        %add one to number of iterations
        numIter = numIter + 1;
        
        %save input PSF sigma in new variable and empty psfSigma for estimation
        psfSigma0 = psfSigma;
        psfSigma = [];
        
        %calculate some numbers that get repeated many times
        psfSigma5 = ceil(5*psfSigma0);
        
        %initialize progress display
        switch numIter
            case 1
                progressText(0,'Estimating PSF sigma');
            otherwise
                progressText(0,'Repeating PSF sigma estimation');
        end
        
        %go over images
        for iImage = goodImages
            
            %read raw image
            imageRaw = movie(:,:,:,iImage);
            
            %get feature positions
            featPos = localMaxima(iImage).cands;
            
            %retain only features that are more than 5*psfSigma0 away from boundaries
            feat2use = find(featPos(:,1) > psfSigma5(1) & featPos(:,1) < imageSizeX - psfSigma5(1) & featPos(:,2) > psfSigma5(1) & featPos(:,2) < imageSizeY - psfSigma5(1) & featPos(:,3) > psfSigma5(2) & featPos(:,3) < imageSizeZ - psfSigma5(2));
            featPos = featPos(feat2use,:);
            
            %if there is more than one feature ...
            if length(feat2use) > 1
                
                % XY distance matrix
                distXY = createDistanceMatrix(featPos(:,1:2),featPos(:,1:2));
                % Z distance matrix
                distZ = abs(createDistanceMatrix(featPos(:,3),featPos(:,3)));
                % find xy and z distance above min distances
                % 10*psf for XY
                distXY = distXY > ceil(10*psfSigma0(1));
                % 10*psf for Z
                distZ = distZ > ceil(10*psfSigma0(2));
                
                % make the diagonals for each matrix true
                for iLoop = 1:size(distXY,1)
                    distXY(iLoop,iLoop) = true;
                    distZ(iLoop,iLoop) = true;
                end
                
                % now find rows that have all true in XY OR Z
                feat2use = all(distXY | distZ, 2);
                
                featPos = featPos(feat2use,:);
            end
            
            %go over the selected features and estimate psfSigma
            numFeats = size(featPos,1);
            
            parameters = zeros(numFeats,7);
            if numFeats >= 1
                
                for iFeat = 1 : numFeats
                    
                    %crop image around selected feature
                    lowerBoundXY = featPos(iFeat,1:2) - psfSigma5(1);
                    upperBoundXY = featPos(iFeat,1:2) + psfSigma5(1);
                    lowerBoundZ = featPos(iFeat,3) - psfSigma5(2);
                    upperBoundZ = featPos(iFeat,3) + psfSigma5(2);
                    imageCropped = imageRaw(lowerBoundXY(1):upperBoundXY(1),lowerBoundXY(2):upperBoundXY(2),lowerBoundZ:upperBoundZ,1);
                    
                    %estimate sigma if image region contains no NaNs
                    %NaNs appear due to cropping
                    if all(~isnan(imageCropped(:)))
                        
                        % initial guess of parameters, leave amp and bg to
                        % be estimated
                        initGuess = [psfSigma5(1)+1 psfSigma5(1)+1 psfSigma5(2)+1 NaN psfSigma0(1) psfSigma0(2) NaN];
                        
                        %fit image and estimate sigma of Gaussian
                        parameters(iFeat,:) = GaussFitND(imageCropped,[],fitParameters,initGuess);
                        
                    else %otherwise assign NaN
                        
                        parameters(iFeat,:) = NaN;
                        
                    end
                    
                end
                
                %add to array of sigmas
                psfSigma = [psfSigma; parameters(:,5:6)]; %#ok<AGROW>
                
            end %(if numFeats >= 1)
            
            %display progress
            switch numIter
                case 1
                    progressText(iImage/max(goodImages),'Estimating PSF sigma');
                otherwise
                    progressText(iImage/max(goodImages),'Repeating PSF sigma estimation');
            end
            
        end %(for iImage = images2use)
        
        %estimate psfSigma as the robust mean of all the sigmas from the fits
        %get rid of NaNs from cropped regions
        if ~isempty(psfSigma)
            psfSigma = psfSigma((~isnan(psfSigma(:,1)) & ~isnan(psfSigma(:,2))),:);
        end
        numCalcs = size(psfSigma,1);
        if numCalcs > 0
            
            [psfSigmaXY_temp,~,inlierIndxXY_temp] = robustMean(psfSigma(:,1));
            [psfSigmaZ_temp,~,inlierIndxZ_temp] = robustMean(psfSigma(:,2));
            
            psfSigma = [psfSigmaXY_temp psfSigmaZ_temp];
            
            acceptCalc = (numCalcs >= 100 && length(inlierIndxXY_temp) >= 0.5*numCalcs && length(inlierIndxZ_temp) >= 0.5*numCalcs) || ...
                (numCalcs >= 50 && length(inlierIndxXY_temp) >= 0.7*numCalcs && length(inlierIndxZ_temp) >= 0.7*numCalcs) || ...
                (numCalcs >= 10 && length(inlierIndxXY_temp) >= 0.9*numCalcs && length(inlierIndxZ_temp) >= 0.9*numCalcs);
            
        else
            
            acceptCalc = 0;
            
        end
        
        %show new sigma if estimation is accepted
        if acceptCalc
            disp(sprintf('PSF sigma = [%1.3f %1.3f] (%d XY inliers and %d Z inliers out of %d observations)',...
                psfSigma(1),psfSigma(2),length(inlierIndxXY_temp),length(inlierIndxZ_temp),numCalcs));
        else %otherwise alert user that input sigma was retained
            psfSigma = psfSigmaIn;
            disp('Not enough observations to change PSF sigma, using input PSF sigma');
        end
        
    end %(while numIter <= numSigmaIter && acceptCalc && ((psfSigma-psfSigma0)/psfSigma0 > 0.05))
    
    %if maximum number of iterations has been performed but sigma value is not converging
    if numIter == numSigmaIter+1 && acceptCalc && any((abs(psfSigma - psfSigma0) ./ psfSigma0) > [0.05 0.05])
        psfSigma = psfSigmaIn;
        disp('Estimation terminated (no convergence), using input PSF sigma');
    end
    
end %(if numSigmaIter)

%% Mixture-model fitting

% save the pixelSize
pixelSize = [dataStruct.dataProperties.PIXELSIZE_XY dataStruct.dataProperties.PIXELSIZE_XY dataStruct.dataProperties.PIXELSIZE_Z];

%initialize initCoord
initCoord(1:nFrames) = struct('allCoord',[],'allCoordPix',[],'nSpots',0,'amp',[],'bg',[]);

progressText(0,'Mixture-model fitting');

%go over all non-empty images ...
for iImage = goodImages
    
    
    %get raw image
    imageRaw = movie(:,:,:,iImage);
    
    % get cands
    cands = localMaxima(iImage).cands;
    
    try %try to detect features in this frame
        
        %fit with mixture-models
        
        % first jointly fit all with no N+1 to remove false spots
        %coordList = spotMMFit(imageRaw,cands,psfSigma,'fitNPlusOne',0,'influenceRadius',1.7*psfSigma,'overlapFactor',psfSigma);
        
        % next N+1 fit with clustering, and don't perform an amp test
        %coordList = spotMMFit(imageRaw,coordList(:,1:3),psfSigma,'fitNPlusOne',1,'influenceRadius',6*psfSigma,'overlapFactor',13*psfSigma,'amplitudeCutoff',1);
        
        % final full joint fit to remove bad spots
        %[coordList,ampList,bgList] = spotMMFit(imageRaw,coordList(:,1:3),psfSigma,'fitNPlusOne',0,'influenceRadius',1.7*psfSigma,'overlapFactor',1000*psfSigma);
        %[coordList,ampList,bgList] = spotMMFit(imageRaw,cands,psfSigma,'fitNPlusOne',1,'influenceRadius',6*psfSigma,'overlapFactor',11*psfSigma);
        
        [coordList,ampList,bgList] = fastGaussFit( imageRaw,cands,psfSigma );
        
        %save results
        initCoord(iImage).nSpots = size(coordList,1);
        initCoord(iImage).allCoordPix = coordList;
        initCoord(iImage).amp = ampList;
        initCoord(iImage).bg = bgList;
        
        % calc real space coordinates
        initCoord(iImage).allCoord = initCoord(iImage).allCoordPix .* repmat(pixelSize,initCoord(iImage).nSpots,2);
        
        %check whether frame is empty
        if initCoord(iImage).nSpots == 0
            emptyFrames = [emptyFrames; iImage]; %#ok<AGROW>
        end
        
    catch %#ok<CTCH> %if detection fails
        
        %label frame as empty
        emptyFrames = [emptyFrames; iImage]; %#ok<AGROW>
        
        %add this frame to the array of frames with failed mixture-model
        %fitting
        framesFailedMMF = [framesFailedMMF; iImage]; %#ok<AGROW>
        
    end
    
    %display progress
    progressText(iImage/length(goodImages),'Mixture-model fitting');
end

%% Post-processing

%sort list of empty frames, keep only unique frames
emptyFrames = unique(emptyFrames);

%store empty frames and frames where detection failed in structure
%exceptions
exceptions = struct('emptyFrames',emptyFrames,'framesFailedLocMax',framesFailedLocMax,'framesFailedMMF',framesFailedMMF');

% save results
initCoord(1).exceptions = exceptions;
initCoord(1).localMaxima = localMaxima;
dataStruct.dataProperties.psfSigma = psfSigma;
dataStruct.initCoord = initCoord;


% initCoord depends on dataProperites
dependencies = struct('dataProperties',[]);
dataPropName = dataStruct.dataPropertiesName;
dataPropV = getVersion(dataPropName);
dependencies.dataProperties = dataPropV;

% save into the first initCoords
dataStruct.initCoord(1).dependencies = dependencies;

%go back to original warnings state
warning(warningState);


%% ~~~ the end ~~~

