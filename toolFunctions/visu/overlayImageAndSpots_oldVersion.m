function overlayImageAndSpots_oldVersion( dataStruct )
%OVERLAYIMAGEANDSPOTS Overlays slices of an image and overlays spot detection results
%   EHarry 2012

global imageAll


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

[~,~,ext] = fileparts(imageInfo.imageName);

switch ext
    case '.mat'
        imageInfo.imageType = 'matFile';
    case '.dv'
        imageInfo.imageType = 'dv';
    otherwise
        error('--overlayImageAndSpots: Error, unreconised movie type');
end

switch imageInfo.imageType
    
    case 'matFile'
        imageAll = readOMEMatFile(imageInfo.imageName,[],1,imageInfo.imageDecon,imageInfo.imageCrop);
        
    case 'dv'
        imageAll = cdLoadMovie({imageInfo.imageName,'raw'},[],struct('crop',imageInfo.imageCrop,'maxSize',1e20));
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
        
        image = imageAll(:,:,:,currentTimePoint);
        
%         switch imageInfo.imageType
%             
%             case 'matFile'
%                 image = readOMEMatFile(imageInfo.imageName,currentTimePoint,1,imageInfo.imageDecon,imageInfo.imageCrop);
%                 
%             case 'dv'
%                 image = cdLoadMovie({imageInfo.imageName,'raw'},[],struct('frames2load',{{currentTimePoint}},'crop',imageInfo.imageCrop,'maxSize',1e20));
%         end
        
        maxIn = max(image(:));
        minIn = min(image(:));
        image = image(:,:,currentZSec);
        spots = imageInfo.spots(currentTimePoint).allCoordPix;
        
        spots = spots((spots(:,3) < (currentZSec+0.5) & spots(:,3) > (currentZSec-0.5)),:);
        
        set(figH,'Name', sprintf('Spots %s (t = %2d / %2d, z = %2d / %2d)',imageInfo.projectName,currentTimePoint, timeInfo.nTimePoints,currentZSec,zInfo.nZSec));
        
        imshow(image,[minIn maxIn]);
        if ~initialPlot
            %set(gca,cameraPropertyNames,cameraProps);
            zoom reset
            set(gca,{'xlim','ylim'},L)
        end
        
        if ~isempty(spots)
            hold on
            plot(spots(:,2),spots(:,1),'g*');
            hold off
        end
        
        % end of function plotPlanes
        
        
    end
end

