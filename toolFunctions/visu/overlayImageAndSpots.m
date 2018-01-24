function overlayImageAndSpots( dataStruct )
%OVERLAYIMAGEANDSPOTS Overlays slices of an image and overlays spot detection results
%   EHarry 2012

global imageAll
global lastFrame
lastFrame = 0;
global lastImage
global lastMaxIn
global lastMinIn
global lastSpotImage
global lastColorImage
global lastMap
global lastSpots
global coordsAll
global spotImageAll
global colorImageAll
global mapAll

% check input for empty
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end
% check whether analysis has been done
if isempty(dataStruct.initCoord)
    return
end

nTimePoints = dataStruct.dataProperties.movieSize(4);
nZSec = dataStruct.dataProperties.movieSize(3);
imageInfo.spots = dataStruct.initCoord;
imageInfo.projectName = dataStruct.projectName;
imageInfo.imageName = fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName);

if isfield(dataStruct.dataProperties,'decon')
    imageInfo.imageDecon = dataStruct.dataProperties.decon;
else
    imageInfo.imageDecon = 0;
end

if isfield(dataStruct.dataProperties,'crop')
    imageInfo.imageCrop = dataStruct.dataProperties.crop;
else
    imageInfo.imageCrop = [];
end

if isfield(dataStruct.dataProperties,'psfSigma')
    imageInfo.psfSigma = dataStruct.dataProperties.psfSigma;
else
    if isfield(dataStruct.dataProperties,'FILTERPRM')
        imageInfo.psfSigma = [dataStruct.dataProperties.FILTERPRM(1) dataStruct.dataProperties.FILTERPRM(3)];
    else
        imageInfo.psfSigma = [1 1];
    end
end

[~,~,ext] = fileparts(imageInfo.imageName);

switch ext
    case '.mat'
        imageInfo.imageType = 'matFile';
    case '.dv'
        imageInfo.imageType = 'dv';
    otherwise
        error('--overlayImageAndSpots: Error, unreconised movie type');
end

% initialize plot window
figH = figure('Name',sprintf('Spots %s (t = %2d / %2d, z = %2d / %2d)',dataStruct.projectName,1, nTimePoints,1,nZSec));
set(figH,'KeyPressFcn',@figure_keyPress);
% user data of the figure stores the current time point and the number of frames
imageInfo.timeInfo = struct('nTimePoints',nTimePoints,...
    'currentTimePoint',1);
imageInfo.zInfo = struct('nZSec',nZSec,'currentZSec',1);
set(figH,'UserData',imageInfo);
plotPlanes(figH,1);


%% LOCAL FUNCTIONS

    function figure_keyPress(src,event)
        
        imageInfo = get(src,'UserData');
        timeInfo = imageInfo.timeInfo;
        zInfo = imageInfo.zInfo;
        
        switch event.Key
            case 'rightarrow'
                timeInfo.currentTimePoint = timeInfo.currentTimePoint + 1;
            case 'leftarrow'
                timeInfo.currentTimePoint = timeInfo.currentTimePoint - 1;
            case 'uparrow'
                zInfo.currentZSec = zInfo.currentZSec + 1;
            case 'downarrow'
                zInfo.currentZSec = zInfo.currentZSec - 1;
        end
        
        if timeInfo.currentTimePoint > timeInfo.nTimePoints
            timeInfo.currentTimePoint = 1;
        end
        if timeInfo.currentTimePoint < 1
            timeInfo.currentTimePoint = timeInfo.nTimePoints;
        end
        if zInfo.currentZSec > zInfo.nZSec
            zInfo.currentZSec = 1;
        end
        if zInfo.currentZSec < 1
            zInfo.currentZSec = zInfo.nZSec;
        end
        
        imageInfo.timeInfo = timeInfo;
        imageInfo.zInfo = zInfo;
        
        set(src,'UserData',imageInfo);
        plotPlanes(src,0,event.Key);
        
    end % end of function figure_keyPress

    function plotPlanes(figH,initialPlot,key)
        
        %cameraPropertyNames = {'CameraPosition'};%, 'CameraTarget', 'CameraUpVector', 'CameraViewAngle'};
        
        
        % figure is still open
        figure(figH)
        %cameraProps = get(gca,cameraPropertyNames);
        L = get(gca,{'xlim','ylim'});
        
        imageInfo = get(figH,'UserData');
        timeInfo = imageInfo.timeInfo;
        currentTimePoint = timeInfo.currentTimePoint;
        zInfo = imageInfo.zInfo;
        currentZSec = zInfo.currentZSec;
        %         planeData = imageInfo.planeFit(currentTimePoint);
        
        if initialPlot
            
            disp('Loading images...');
            switch imageInfo.imageType
                
                case 'matFile'
                    imageAll = readOMEMatFile(imageInfo.imageName,[],1,imageInfo.imageDecon,imageInfo.imageCrop);
                    
                case 'dv'
                    imageAll = cdLoadMovie({imageInfo.imageName,'raw'},[],struct('crop',imageInfo.imageCrop,'maxSize',1e20));
            end
            
            [sizeX,sizeY,sizeZ,~] = size(imageAll);
            
            %             coordsAll(1:timeInfo.nTimePoints) = struct('coords',[]);
            %             spotImageAll(1:timeInfo.nTimePoints) = struct('spotImage',[]);
            %             colorImageAll(1:timeInfo.nTimePoints) = struct('colorImage',[]);
            %             mapAll(1:timeInfo.nTimePoints) = struct('map',[]);
            %
            progressText(0,'Calculating spot pixel intnsities...');
            for t = 1:timeInfo.nTimePoints
                
                coords = [imageInfo.spots(t).allCoordPix(:,1:3),imageInfo.spots(t).amp(:,1)];
                %                 [ spotImage, colorImage, map ] = makeColorOverlayedImage( coords, imageInfo.psfSigma, sizeX , sizeY, sizeZ );
                [spotImage, colorImage, map] = makeColorOverlayedImage( coords, imageInfo.psfSigma, sizeX , sizeY, sizeZ );
                
                coordsAll(t).coords = coords;
                spotImageAll(t).spotImage = spotImage;
                colorImageAll(t).colorImage = colorImage;
                mapAll(t).map = map;
                
                progressText(t/timeInfo.nTimePoints,['Calculating spot pixel intnsities at t=',int2str(t)]);
            end
        end
        
        
        
      %  if lastFrame ~= currentTimePoint
            
            
            %             switch imageInfo.imageType
            %
            %                 case 'matFile'
            %                     image = readOMEMatFile(imageInfo.imageName,currentTimePoint,1,imageInfo.imageDecon,imageInfo.imageCrop);
            %
            %                 case 'dv'
            %                     image = cdLoadMovie({imageInfo.imageName,'raw'},[],struct('frames2load',{{currentTimePoint}},'crop',imageInfo.imageCrop,'maxSize',1e20));
            %             end
            
            image = imageAll(:,:,:,currentTimePoint);
            
            
            maxIn = max(image(:));
            minIn = min(image(:));
            %         image = image(:,:,currentZSec);
            %         spots = imageInfo.spots(currentTimePoint).allCoordPix;
            %             coords = [imageInfo.spots(currentTimePoint).allCoordPix(:,1:3),imageInfo.spots(currentTimePoint).amp(:,1)];
            %             [ spotImage, colorImage, map ] = makeColorOverlayedImage( coords, imageInfo.psfSigma, size(image,1) , size(image,2), size(image,3) );
            
            coords = coordsAll(currentTimePoint).coords;
            spotImage = spotImageAll(currentTimePoint).spotImage;
            colorImage = colorImageAll(currentTimePoint).colorImage;
            map = mapAll(currentTimePoint).map;
            
            lastFrame = currentTimePoint;
            lastImage = image;
            lastMaxIn = maxIn;
            lastMinIn = minIn;
            lastSpotImage = spotImage;
            lastColorImage = colorImage;
            lastMap = map;
            lastSpots = coords;
            
%         else
%             
%             image = lastImage;
%             maxIn = lastMaxIn;
%             minIn = lastMinIn;
%             spotImage = lastSpotImage;
%             colorImage = lastSpotImage;
%             map = lastMap;
%             coords = lastSpots;
%             
%         end
%         
        
        
        
        %         spots = spots((spots(:,3) < (currentZSec+0.5) & spots(:,3) > (currentZSec-0.5)),:);
        index = find(coords(:,3) < (currentZSec+0.5) & coords(:,3) > (currentZSec-0.5));
        
        set(figH,'Name', sprintf('Spots %s (t = %2d / %2d, z = %2d / %2d)',imageInfo.projectName,currentTimePoint, timeInfo.nTimePoints,currentZSec,zInfo.nZSec));
        
        %         imshow(image,[minIn maxIn]);
        imshow(squeeze(max(image(:,:,currentZSec),[],3)),[minIn maxIn])
        hold on
        h = imshow(squeeze(colorImage(:,:,currentZSec,:)));
        temp = squeeze(spotImage(:,:,currentZSec));
        temp = temp - min(temp(:));
        temp = temp ./ max(temp(:));
        set(h, 'AlphaData', temp);
        hold off
        
        if ~initialPlot
            %set(gca,cameraPropertyNames,cameraProps);
            zoom reset
            set(gca,{'xlim','ylim'},L)
        end
        
        if ~isempty(index)
            hold on
            %             plot(spots(:,2),spots(:,1),'g*');
            
            % add text
            for i = index'
                
                textCoord = coords(i,1:2) + [2 -2];
                col = map(i,:);
                text(textCoord(2),textCoord(1),['\color[rgb]{',num2str(col(1)),' ',num2str(col(2)),' ',num2str(col(3)),'}',int2str(i)])
                plot(coords(i,2),coords(i,1),'Color',col,'Marker','*')
                
            end
            
            hold off
        end
        
        % end of function plotPlanes
        
        
    end
end

