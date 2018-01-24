function dataStruct = makiInitCoord2D(dataStruct, verbose, movieType, cutoffFix)

% EHarry October 2011

%% ORIGIANL HEADER
% % %MAKIINITCOORD finds initial guesses for mammalian kinetochores
% % %
% % % SYNOPSIS: dataStruct = makiInitCoord(dataStruct)
% % %
% % % INPUT dataStruct: dataStruct as in makimakeDataStruct with at least the
% % %                       fields
% % %                         .rawMovieName - may contain the actual movie
% % %                         .rawMoviePath
% % %                         .dataProperties
% % %                       Alternatively, a makiData object can be passed.
% % %                   verbose (opt): 0: plot nothing
% % %                                  1: show progressText (default)
% % %                                  2: + show cutoff-plots
% % %                                  3: + show initial coordinate guesses
% % %                                     ordered by amplitude for every frame
% % %                   movieType: 1 for DV movies, 2 for STK movies.
% % %                              Optional. Default: 1.
% % %                   cutoffFix: 2-element vector indicating
% % %                           cutoff-type/cutoff-value for the selection of
% % %                           the "good" local maxima. Default: []
% % %                           This input is to be used for testing purposes
% % %                           only; a pre-set cutoff value should be stored
% % %                           in dataStruct.dataProperties
% % %
% % % OUTPUT dataStruct.initCoords: Structure of length nTimepoints with fields
% % %                       .allCoord     [x y z sx sy sz] coords and sigmas in
% % %                           microns, corrected for refocussing. Sigma is
% % %                           0.25 pixels
% % %                       .allCoordPix  same as coords, but in pixels and
% % %                           without correction.
% % %                       .correctionMu z-correction in microns
% % %                       .nSpots       number of spots
% % %                       .initAmp      maximum pixel intensity and estimated
% % %                           local noise of local maxima
% % %                           This tends to be better than .amp
% % %                       .amp          amplitude estimate from local
% % %                           integral around locMax
% % %                       .data4MMF     [x y z] pixel coordinates of two
% % %                           spots right above and right below the cutoff,
% % %                           which can be used with a function like
% % %                           detectSpots_MMF_findAmplitudeCutoff.m to
% % %                           calculate the t-test cutoff for mixture model
% % %                           fitting (see code snipped at end of fcn)
% % %
% % % REMARKS
% % %
% % % created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
% % %
% % % created by: jdorn
% % % DATE: 28-Jun-2007
% % %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(verbose)
    verbose = 1;
end

if nargin < 3 || isempty(movieType)
    movieType = 1;
end

if nargin < 4 || isempty(cutoffFix)
    cutoffFix = [];
end

%=========================
%% COLLECT INPUT & SETUP
%=========================

if isobject(dataStruct)
    % if dataStruct is an object, movie is in imageData
    dataObject = true;
    % change save-mode of initCoord so that it doesn't autosave 100 times
    dataStruct.propertyNames{3,2}.saveMode = 0;
else
    dataObject = false;
    if ischar(dataStruct.rawMovieName)
        rawMovie = fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName);
    else
        rawMovie = [];
    end
end
dataProperties = dataStruct.dataProperties;

% find more parameters
% betterBackground: estimate background after masking signal
if ~isfield(dataProperties,'betterBackground')
    if dataObject
        %betterBackground = false;
        betterBackground = true;
    else
        betterBackground = false;
    end
end
% isNanMask: switch that checks whether there needs to be a correction for
% nanMasks
if ~isfield(dataProperties,'isNanMask')
    if dataObject
        isNanMask = 2;
    else
        isNanMask = 0;
    end
end

% setup filter parameters. Increase speed by doing incremental filtering
% for background
backgroundFilterParms = dataProperties.FT_SIGMA * 14; % 14 is 15-1
backgroundFilterParms(4:6) = roundOddOrEven(backgroundFilterParms,'odd','inf');
% create background mask already
backgroundFilter = GaussMask3D(backgroundFilterParms(1:3),...
    backgroundFilterParms(4:6),[],1,[],[],1);
signalFilter = GaussMask3D(dataProperties.FILTERPRM(1:3),...
    dataProperties.FILTERPRM(4:6),[],1,[],[],1);

% make separated noise mask
noiseMask = {...
    ones(dataProperties.FILTERPRM(4),1,1)./dataProperties.FILTERPRM(4),...
    ones(1,dataProperties.FILTERPRM(5),1)./dataProperties.FILTERPRM(5),...
    ones(1,1,dataProperties.FILTERPRM(6))./dataProperties.FILTERPRM(6),...
    };

% for conversion to microns in the end: get pixelSize
pixelSize = [dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_Z];

% get halfPatchSize to adjust centroid result. The center pixel of a 5x5x5
% array is (3,3,3), thus, we have to subtract that from centroid coordinates
halfPatchSize = dataProperties.FILTERPRM(4:6)/2+0.5;

% setup lists
nTimepoints = dataProperties.movieSize(4);
initCoordRaw = cell(nTimepoints,1); % raw initial coords (before testing)
% set up initCoord-Structure
tmp(1:nTimepoints,1) = struct('allCoord',[],...
    'allCoordPix',[],'correctionMu',[],'nSpots',[],...
    'initAmp',[],'amp',[],'data4MMF',[]);
dataStruct.initCoord = tmp;

% missing frames are indicated by empty cropMasks. Don't analyze them!
% For safety reason, this check is only performed on data objects
goodTimes = 1:nTimepoints;
if dataObject
    % since this is dataObj only, no need to worry about compatibility
    if dataProperties.initCoord.rmEmptyMasks
        goodTimes = find(cellfun(@(x)(~isempty(x)),...
            dataStruct.imageData.cropInfo.cropMask))';
    end
end

% read minimum number of requested spots per frame. If not set, assume 20
if dataObject
    if isfield(dataProperties.initCoord,'minSpotsPerFrame')
        minSpotsPerFrame = dataProperties.initCoord.minSpotsPerFrame;
    else
        minSpotsPerFrame = 20;
    end
else
    minSpotsPerFrame = 20;
end

 blurKernelHigh = fspecial('gaussian', 21, 4);
 blurKernelLow  = fspecial('gaussian', 21, 1);


%=====================
%% MAIN LOOP
%=====================
if verbose
    progressText;
end
nTimes = max(goodTimes); % not the best way to count, but it works.
for t=goodTimes
    
    % -- load images --
    % raw = current raw movie frame
    % filtered = current PSF-filtered movie frame
    % then calculate
    % background = current 10xPSF-filtered movie frame
    % amplitude = filtered - background
    % noise = locAvg(var(raw-filtered))
    
    if dataObject
        % isNanMask == 2 indicates that we should crop only later in order
        % to increase processing speed.
        raw = dataStruct.imageData.getFrame(t,[],[],isNanMask==1);
    elseif isempty(rawMovie)
        % movie has been passed directly
        raw = dataStruct.rawMovieName(:,:,:,t); %EH changed from 5d to 4d 09/06/2011
        cropInfo = dataProperties.crop;
        if ~isempty(cropInfo)
            cropInfo = cropInfo(:,1:2);
            cropInfo(1,:) = max(cropInfo(1,:),[1 1]); %defaults for min
            if cropInfo(2,1) == 0 %defaults for max
                cropInfo(2,1) = size(raw,1);
            end
            if cropInfo(2,2) == 0
                cropInfo(2,2) = size(raw,2);
            end
            raw = raw(cropInfo(1,1):cropInfo(2,1),cropInfo(1,2):cropInfo(2,2),:); %crop image
        end
    else
        switch movieType
            case 1 %DV files
                raw = cdLoadMovie({rawMovie,'raw'},[],struct('frames2load',{{t}},...
                    'crop',dataProperties.crop,'maxSize',dataProperties.maxSize));
            case 2 %STK files
                movieTPName = [rawMovie(1:end-5) num2str(t) '.STK']; %name of files storing current time point
                imageStack = metaTiffRead(movieTPName,[],[],0); %read stack using Nedelec's code
                raw = double(cat(3,imageStack.data)); %extract image information
                cropInfo = dataProperties.crop; %cropping information
                if ~isempty(cropInfo)
                    cropInfo = cropInfo(:,1:2);
                    cropInfo(1,:) = max(cropInfo(1,:),[1 1]); %defaults for min
                    if cropInfo(2,1) == 0 %defaults for max
                        cropInfo(2,1) = size(raw,1);
                    end
                    if cropInfo(2,2) == 0
                        cropInfo(2,2) = size(raw,2);
                    end
                    raw = raw(cropInfo(1,1):cropInfo(2,1),cropInfo(1,2):cropInfo(2,2),:); %crop image
                end
        end
    end
    %remove offset
    offset = min(raw(:));
    raw = raw -  offset;
    % filter movie with gauss-filter
    
    if betterBackground
        % there will be NaNs in the masked image. Therefore, request
        % filtered data from imageDataObject. Ideally, loadType is set to
        % fcn or reqKeep
        
        filtered = fastGauss3D(raw,[],dataProperties.FILTERPRM(4:6),2-(isNanMask==1),signalFilter,1);
        % filter signal first, then run smaller filter with background to
        % save some time
        
        background = fastGauss3D(filtered,[],backgroundFilterParms(4:6),2-(isNanMask==1),backgroundFilter,1);
        
        % mask signal pixels, then recalculate background
        sigMask = raw>background;
        % remove some of the spurious hits. Use 2d mask for speed and
        % memory (3d strel won't work on binary image with imopen)
        sigMask = imopen(sigMask,strel('disk',1));
        rawMsk = raw.*~sigMask;
        clear sigMask
        % fill holes. ConvNan can in principle handle holes, but they might
        % be on the large side. Similarly, if fillZeroHoles misses a few
        % zeros, convNan can help out
        for z=1:dataProperties.movieSize(3)
            rawMsk(:,:,z)=blkproc(rawMsk(:,:,z),[21 21],@fillZeroHoles);
        end
        rawMsk(rawMsk==0) = NaN;
        background = fastGauss3D(rawMsk,[],dataProperties.FILTERPRM(4:6),1,signalFilter,1);
        %background = convNan(rawMsk,backgroundFilter,backgroundFilterParms(4:6),1);
        
    else
        % only use NaN-image-reduction and new addBorders for dataObj to guarantee backward
        % compatibility
        %filtered = fastGauss3D(raw,[],dataProperties.FILTERPRM(4:6),2-(isNanMask==1),signalFilter,dataObject);
        % filtering with sigma=15 is equal to filtering with 1 and then
        % with 14
        %background = fastGauss3D(filtered,[],backgroundFilterParms(4:6),2-(isNanMask==1),backgroundFilter,dataObject);
    end
    
    % amplitude is filtered image - background. This underestimates the
    % true signal. The underestimation becomes stronger if betterBackground
    % is not used.
    %amplitude = filtered - background;
    
    startFrame = 1;
    endFrame = size(raw,3);
    
    [movieInfo(1:endFrame,1).xCoord] = deal([]);
    [movieInfo(1:endFrame,1).yCoord] = deal([]);
    [movieInfo(1:endFrame,1).amp] = deal([]);
    [movieInfo(1:endFrame,1).int] = deal([]);
    
    noise = raw;
    
    for iFrame = startFrame:endFrame
        
        
        img = raw(:,:,iFrame);
        [imL,imW] = size(img);
        lowPass = imfilter(img,blurKernelLow);
        highPass = imfilter(img,blurKernelHigh);
        filterDiff = lowPass-highPass;
        filterDiff(1:10,:)=0;
        filterDiff(end-10:end,:)=0;
        filterDiff(:,1:10)=0;
        filterDiff(:,end-10:end)=0;
        % thickness of intensity slices is average std from filterDiffs over
        % from one frame before to one frame after
        stepSize=std(filterDiff(:));
        thresh=3*stepSize;
        % we assume each step size down the intensity profile should be on
        % the order of the size of the background std; here we find how many
        % steps we need and what their spacing should be. we also assume peaks
        % should be taller than 3*std
        
        % noise is local average (averaged over filter support) of squared
        % residuals of raw-filtered image
        noise(:,:,iFrame) = (raw(:,:,iFrame)-filterDiff).^2;
       % noise(:,:,iFrame) = fastGauss3D(noise(:,:,iFrame),[],dataProperties.FILTERPRM(4:5),2-(isNanMask==1),noiseMask);
        noise(:,:,iFrame) = imfilter(noise(:,:,iFrame),blurKernelLow);
        
        nSteps = round((nanmax(filterDiff(:))-thresh)/(stepSize));
        threshList = linspace(nanmax(filterDiff(:)),thresh,nSteps);
        % compare features in z-slices startest from the highest one
        slice2 = zeros(imL,imW);
        for p = 1:length(threshList)-1
            % slice1 is top slice; slice2 is next slice down
            % here we generate BW masks of slices
            if p==1
                slice1 = filterDiff>threshList(p);
            else
                slice1 = slice2;
            end
            slice2 = filterDiff>threshList(p+1);
            % now we label them
            featMap1 = bwlabel(slice1);
            featMap2 = bwlabel(slice2);
            featProp2 = regionprops(featMap2,'PixelIdxList');
            % loop thru slice2 features and replace them if there are 2 or
            % more features from slice1 that contribute
            for iFeat = 1:max(featMap2(:))
                pixIdx = featProp2(iFeat,1).PixelIdxList; % pixel indices from slice2
                featIdx = unique(featMap1(pixIdx)); % feature indices from slice1 using same pixels
                featIdx(featIdx==0) = []; % 0's shouldn't count since not feature
                if length(featIdx)>1 % if two or more features contribute...
                    slice2(pixIdx) = slice1(pixIdx); % replace slice2 pixels with slice1 values
                end
            end
        end
        % label slice2 again and get region properties
        featMap2 = bwlabel(slice2);
        featProp2 = regionprops(featMap2,'PixelIdxList','Area');
        % here we sort through features and retain only the "good" ones
        % we assume the good features have area > 2 pixels
        goodFeatIdx = find(vertcat(featProp2(:,1).Area)>2);
        %    goodFeatIdxI = find(vertcat(featProp2(:,1).MaxIntensity)>2*cutOffValueInitInt);
        %    goodFeatIdx = intersect(goodFeatIdxA,goodFeatIdxI);
        % make new label matrix and get props
        featureMap = zeros(imL,imW);
        featureMap(vertcat(featProp2(goodFeatIdx,1).PixelIdxList)) = 1;
        [featMapFinal,nFeats] = bwlabel(featureMap);
        verDate=version('-date');
        if str2double(verDate(end-3:end))>=2008
            featPropFinal = regionprops(featMapFinal,filterDiff,'PixelIdxList','Area','WeightedCentroid','MaxIntensity'); %'Extrema'
        else
            featPropFinal = regionprops(featMapFinal,'PixelIdxList','Area','Centroid');
            for iFeat=1:length(featPropFinal)
                featPropFinal(iFeat,1).WeightedCentroid=featPropFinal(iFeat,1).Centroid; % centroid's close enough...
                featPropFinal(iFeat,1).MaxIntensity=max(filterDiff(featPropFinal(iFeat,1).PixelIdxList)); % find maximum intensity
            end
        end
        if nFeats==0
            yCoord = [];
            xCoord = [];
            amp = [];
            featI = [];
        else
            % centroid coordinates with 0.5 uncertainties for Khuloud's tracker
            yCoord = 0.5*ones(nFeats,2);
            xCoord = 0.5*ones(nFeats,2);
            temp = vertcat(featPropFinal.WeightedCentroid);
            yCoord(:,1) = temp(:,1);
            xCoord(:,1) = temp(:,2);
            % area
            featArea = vertcat(featPropFinal(:,1).Area);
            amp = zeros(nFeats,2);
            amp(:,1) = featArea;
            % intensity
            featInt = vertcat(featPropFinal(:,1).MaxIntensity);
            featI = zeros(nFeats,2);
            featI(:,1) = featInt;
        end
        % make structure compatible with Khuloud's tracker
        movieInfo(iFrame,1).xCoord = xCoord;
        movieInfo(iFrame,1).yCoord = yCoord;
        movieInfo(iFrame,1).amp = amp;
        movieInfo(iFrame,1).int = featI;
       
       
   
    end
    
    
    if isNanMask == 2
        % apply nan-mask now
        if ~isempty(dataStruct.imageData.cropInfo)
            if iscell(dataStruct.imageData.cropInfo.cropMask)
                msk = ~dataStruct.imageData.cropInfo.cropMask{t};
            else
                msk = ~dataStruct.imageData.cropInfo.cropMask;
            end
            msk = repmat(msk,[1,1,dataProperties.movieSize(3)]);
            filtered(msk) = NaN;
        end
    end
    
    %     figure,imshow(max(amplitude,[],3)-mean(amplitude,3),[]);
    %     error('bam')
    
    % find local maxima
%     locMax = loc_max3Df(filtered,[3 3 3]);
%     locMaxIdx = sub2ind(dataProperties.movieSize(1:3),...
%         locMax(:,1),locMax(:,2),locMax(:,3));
    
    % read amplitude and noise at local maxima
%     initCoordTmp = [locMax, amplitude(locMaxIdx), noise(locMaxIdx), ...
%         zeros(length(locMaxIdx),1)];

initCoordTmp = [];


for i =1:length(movieInfo)
    
    if ~isempty(movieInfo(i).xCoord)
        nseIdx = [ceil(movieInfo(i).xCoord(:,1)),ceil(movieInfo(i).yCoord(:,1)),i.*ones(size(movieInfo(i).xCoord,1),1)];
        
        noiseTemp = zeros(size(nseIdx,1),1);
        for j =1:size(nseIdx,1)
            noiseTemp(j) = noise(nseIdx(j,1),nseIdx(j,2),nseIdx(j,3));
        end
        
        initCoordTmp = [initCoordTmp;movieInfo(i).xCoord(:,1),movieInfo(i).yCoord(:,1),ones(size(movieInfo(i).xCoord,1),1).*i,movieInfo(i).int(:,1),noiseTemp, ...
            zeros(size(movieInfo(i).xCoord,1),1)];
    end
end
    initCoordTmp = sortrows(initCoordTmp,-4);
    % only take MAXSPOTS highest amplitudes. This will leave plenty of noise
    % spots, but it will avoid huge arrays
    initCoordTmp = initCoordTmp(1:min(dataProperties.MAXSPOTS,size(initCoordTmp,1)),:);
    
    %remove any spots with negative amplitude, because those are false
    %positives for sure - KJ
    if movieType == 2
        initCoordTmp = initCoordTmp(initCoordTmp(:,4)>0,:);
    end
    
    % loop through all to get sub-pixel positions.
%     raw = raw - background; %overwrite raw to save memory
%     for iSpot = 1:size(initCoordTmp,1)
%         
%         % read volume around coordinate
%         patch = stamp3d(...
%             raw,dataProperties.FILTERPRM(4:6),...
%             initCoordTmp(iSpot,1:3),0);
%         
%         %HLE,KJ - calculate low-index edge patch adjustment if relevant
%         edgeAdjustTmp = initCoordTmp(iSpot,1:3) - halfPatchSize;
%         edgeAdjustTmp = abs(edgeAdjustTmp) .* edgeAdjustTmp<0;
%         
%         % subpixel coord is integer coord plus centroid (subtract
%         % halfPatchSize so that center coordinate of patch is (0,0,0))
%         initCoordTmp(iSpot,1:3) = ...
%             initCoordTmp(iSpot,1:3) + ...
%             centroid3D(patch) - halfPatchSize;
%         
%         %HLE,KJ - correct position due to patch truncation due to proximity
%         %to a low-index edge
%         initCoordTmp(iSpot,1:3) = initCoordTmp(iSpot,1:3) + edgeAdjustTmp;
%         
%         % %         %KJ - Hunter's original attempt - has a bug:
%         % %         %last line must have floor(edgeAdjustTmp)
%         % %         %HLE - Adjust for edge effect. When patch is truncated by low-index
%         % %         %edge of image, the halfPatchSize needs to be adjusted
%         % %         edgeAdjustTmp = initCoordTmp(iSpot,1:3) - ...
%         % %             (dataProperties.FILTERPRM(4:6)-1)./2 ;
%         % %         initCoordTmp(iSpot,1:3) = initCoordTmp(iSpot,1:3) - ...
%         % %             edgeAdjustTmp .* (edgeAdjustTmp < 0);
%         
%         % amplitude guess is integral.
%         initCoordTmp(iSpot,6) = nanmean(patch(:));
%         
%     end
    
    % use maxPix-amplitude to calculate cutoff - meanInt is not consistent
    % with the rest of the measures!
    initCoordTmp = [initCoordTmp,...
        initCoordTmp(:,4)./sqrt(initCoordTmp(:,5)),...
        initCoordTmp(:,4)./sqrt(initCoordTmp(:,5)./max(initCoordTmp(:,4),eps))];
    
    initCoordRaw{t} = initCoordTmp;
    
    if verbose
        progressText(t/nTimes);
    end
    
    if verbose > 2
        % plot frame with initial guesses for coords
        figure('Name',sprintf('frame %i',t))
        imshow(max(amplitude,[],3),[])
        hold on,
        for i=1:size(initCoordTmp,1),
            text(initCoordTmp(i,2),initCoordTmp(i,1),num2str(i),'Color','r');
        end
    end
end % loop timepoints

clear initCoordTmp
allCoord = cat(1,initCoordRaw{:});

%% Find cutoff

% find cutoff based on amplitude/sqrt(noise/amp), though the others are
% very similar. Allow fallback if less than 25 spots per frame are found
% (this indicates that cutFirstHistmode of splitModes failed)
cutoff = zeros(3,1);
cutoff1 = splitModes(allCoord(:,4)); % amplitude
cutoff2 = splitModes(allCoord(:,7)); % amplitude/sqrt(nse) - dark noise
cutoff3 = splitModes(allCoord(:,8)); % amplitude/sqrt(nse/amp) - poisson
if ~isempty(cutoff1)
    cutoff(1) = cutoff1;
else
    cutoff(1) = NaN;
end
if ~isempty(cutoff2)
    cutoff(2) = cutoff2;
else
    cutoff(2) = NaN;
end
if ~isempty(cutoff3)
    cutoff(3) = cutoff3;
else
    cutoff(3) = NaN;
end

% plot all before selecting cutoff so that we can see what went wrong
if verbose > 1
    if dataObject
        figure('Name',sprintf('cutoffs for %s',dataStruct.identifier))
    else
        figure('Name',sprintf('cutoffs for %s',dataStruct.projectName))
    end
    ah(1) = subplot(3,1,1);
    set(ah(1),'NextPlot','add')
    plot(ah(1),[1,nTimepoints],[cutoff(1) cutoff(1)]);
    xlabel('timepoints')
    ylabel('amplitude - signal in grey values')
    ah(2) = subplot(3,1,2);
    set(ah(2),'NextPlot','add')
    plot(ah(2),[1,nTimepoints],[cutoff(2) cutoff(2)]);
    xlabel('timepoints')
    ylabel('amplitude/sqrt(nse) - SNR for dark noise')
    ah(3) = subplot(3,1,3);
    set(ah(3),'NextPlot','add')
    plot(ah(3),[1,nTimepoints],[cutoff(3) cutoff(3)]);
    xlabel('timepoints')
    ylabel('amplitude/sqrt(nse/amp) - SNR for poisson noise')
    for t=goodTimes
        plot(ah(1),t,initCoordRaw{t}(:,4),'+')
        plot(ah(2),t,initCoordRaw{t}(:,7),'+')
        plot(ah(3),t,initCoordRaw{t}(:,8),'+')
    end
end

nn = nan(3,1);
% check for predetermined cutoff
if ~isempty(cutoffFix)
    % store old value
    oldCutoff = cutoff;
    % set cutoffIdx, cutoffCol
    cutoffIdx = cutoffFix(1);
    cutoffCol = cutoffIdx + 5;
    % update cutoff
    cutoff(cutoffIdx) = cutoffFix(2);
else
    
    
    switch movieType
        case 1
            % note: ask only for spots in good frames
            minGood = minSpotsPerFrame*length(goodTimes);
            nn(3) = sum(allCoord(:,8)>cutoff(3))/minGood;
            nn(2) = sum(allCoord(:,7)>cutoff(2))/minGood;
            nn(1) = sum(allCoord(:,4)>cutoff(1))/minGood;
            if nn(3) > 1
                cutoffIdx = 3;
                cutoffCol = 8;
            elseif nn(2) > 1
                cutoffIdx = 2;
                cutoffCol = 7;
            elseif nn(1) > 1
                cutoffIdx = 1;
                cutoffCol = 4;
            else
                error('less than %i spots per frame found. makiInitCoord failed',minSpotsPerFrame)
            end
            
        case 2
            
            %for the HMS data, the signal is very dim and photobleaches a lot,
            %and determining the cutoff on the fly does not work.
            %Thus, I'm fixing the cutoff criterion to the amplitude-to-noise
            %ratio and I'm fixing the cutoff to 2.32 (approximating a 0.01
            %alpha-value in a hypothesis test) - KJ
            cutoff(2) = 2.32;
            cutoffIdx = 2;
            cutoffCol = 7;
            
    end
end

% remember the cutoff criterion used
if ~dataObject
    dataStruct.statusHelp{3,3} = [cutoffIdx,cutoffCol];
end

if verbose > 1
    title(ah(cutoffIdx),'cutoff selected')
    % if we want to look at this, we should do scatterCloud!!
    %     figure('Name',sprintf('cutoff-comparison for %s',dataStruct.projectName))
    %     minC = min(allCoord,[],1);
    %     maxC = max(allCoord,[],1);
    %     subplot(2,2,1)
    %     plot(allCoord(:,4),allCoord(:,6),'.')
    %     hold on, plot([minC(4),maxC(4)],[cutoff2,cutoff2])
    %     hold on, plot([cutoff1, cutoff1],[minC(6),maxC(6)])
    %     subplot(2,2,2)
    %     plot(allCoord(:,7),allCoord(:,6),'.')
    %     hold on, plot([minC(7),maxC(7)],[cutoff2,cutoff2])
    %     hold on, plot([cutoff3, cutoff3],[minC(6),maxC(6)])
    %     subplot(2,2,3)
    %     plot(allCoord(:,4),allCoord(:,7),'.')
    %     hold on, plot([minC(4),maxC(4)],[cutoff3,cutoff3])
    %     hold on, plot([cutoff1, cutoff1],[minC(7),maxC(7)])
end



% loop and store only good locMax. Before that, get z-correction
% get correction from .log file
% load $MOVIENAME.log, parse for start coordinate and determine focus
% adjustment. This should give a z0-value for every frame that we can
% use to correct the um-coords for tracking



for t=goodTimes
    goodIdxL = initCoordRaw{t}(:,cutoffCol) > cutoff(cutoffIdx);
    
    % count good spots
    dataStruct.initCoord(t).nSpots = sum(goodIdxL);
    
    
    
    % store pixel coords. Uncertainty is 0.25 pix
    dataStruct.initCoord(t).allCoordPix = ...
        [initCoordRaw{t}(goodIdxL,1:3),...
        0.25*ones(dataStruct.initCoord(t).nSpots,3)];
    
    % store estimated amplitude and noise
    dataStruct.initCoord(t).initAmp = initCoordRaw{t}(goodIdxL,4:5);
    
    % store integral amplitude
    dataStruct.initCoord(t).amp = [initCoordRaw{t}(goodIdxL,6),...
        zeros(dataStruct.initCoord(t).nSpots,1)];
    
    % store correction
    dataStruct.initCoord(t).correctionMu = 0; % for now
    
    % store coords in microns and correct
    dataStruct.initCoord(t).allCoord = ...
        dataStruct.initCoord(t).allCoordPix.*...
        repmat(pixelSize,dataStruct.initCoord(t).nSpots,2);
    dataStruct.initCoord(t).allCoord(:,3) = ...
        dataStruct.initCoord(t).allCoord(:,3) +...
        dataStruct.initCoord(t).correctionMu;
    
    
    
    
    % store 2 spots above and below cutoff in case we want to get amplitude
    % cutoff for detector later. Note: the spots may be not exactly the two
    % above or below, since we sorted the data according to amplitudes
    % above.
    twoAbove = find(initCoordRaw{t}(:,cutoffCol) > cutoff(cutoffIdx),2,'last');
    twoBelow = find(initCoordRaw{t}(:,cutoffCol) < cutoff(cutoffIdx),2,'first');
    dataStruct.initCoord(t).data4MMF = ...
        initCoordRaw{t}([twoAbove;twoBelow],1:3);
end

% save cutoff info
if dataObject
    if isempty(cutoffFix)
        dataStruct.initCoord(1).cutoff.co = cutoff;
    else
        % with fixed cutoff: store old value
        dataStruct.initCoord(1).cutoff.co = oldCutoff;
    end
    dataStruct.initCoord(1).cutoff.sel = [cutoffIdx,cutoffCol];
    dataStruct.initCoord(1).cutoff.n = nn;
    
    % change save-mode again and save
    dataStruct.propertyNames{3,2}.saveMode = 2;
    dataStruct.initCoord = dataStruct.initCoord;
    
end

% % code to determine optimal amplitude cutoff
% % check if we should just take average testStatistic instead of doing
% % complicated histogram stuff
% cordStruct = makiCoord2Cord(dataStruct.initCoord4MMF);
% [dataProperties] = ...
%     detectSpots_MMF_findAmplitudeCutoff(...
%     rawMovie, cordStruct, dataProperties,[],0)
