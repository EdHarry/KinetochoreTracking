function overlaySistersOnImage( dataStruct, sis2Plot )
%OVERLAYTRACKSONIMAGE Overlay trajData on an image
%   EHarry May 2012

% check input for empty
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end
% check whether analysis has been done
if isempty(dataStruct.initCoord)
    return
end

frames = dataStruct.dataProperties.movieSize(4);
global toPlot
toPlot = round(frames/20);
% if mod(toPlot,2) == 0
%     toPlot = toPlot + 1;
% end
% toPlot = round(frames/2);

if nargin < 2 || isempty(sis2Plot)
    sis2Plot = 1:length(dataStruct.sisterList);
end

global sis2PlotTmp
sis2PlotTmp = sis2Plot;

% load tracks and coords
tracks = dataStruct.tracks;
initCoord = dataStruct.initCoord;

% get sisterTracks
sisTracks = dataStruct.sisterList(1).trackPairs(sis2Plot,1:2);

% get image coordinate of tracks
coords(1:2*size(sisTracks,1)) = struct('imageCoords',[]);
counter = 0;
for i = 1:size(sisTracks,1)
    for ii = 1:2
        imageCoords = NaN(frames,3);
        featIdx = getFeatIdx(tracks(sisTracks(i,ii)),frames);
        for t = 1:frames
            if ~isnan(featIdx(t))
                imageCoords(t,:) = initCoord(t).allCoordPix(featIdx(t),1:3);
            end
        end
        counter = counter + 1;
        coords(counter).imageCoords = imageCoords;
    end
end


global imageAll


nTimePoints = dataStruct.dataProperties.movieSize(4);
nZSec = dataStruct.dataProperties.movieSize(3);


imageInfo.spots = coords;


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
        % make invisibe while processing
        %         set(figH,'Visible','off');
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
        
        
        %spots = imageInfo.spots(currentTimePoint).allCoordPix;
        
        trackCoord = imageInfo.spots;
        
        spots = [];
        trackSeg = [];
        c=0;
        for iTr = 1:length(trackCoord)
            clear trackSegTmp
            c=c+1;
            newSpots = trackCoord(iTr).imageCoords(currentTimePoint,:);
            % if ~isnan(newSpots(1))
            spots = [spots; newSpots];
            
            
            
            numTimePoints = toPlot*2;
            co_minus = -1;
            %co_plus = numTimePoints + 1;
            for j = 1:toPlot
                co_minus = co_minus + 1;
                %co_plus = co_plus - 1;
                
                currentSpotsX = newSpots(1);
                currentSpotsY = newSpots(2);
                currentSpotsZ = newSpots(3);
                if currentTimePoint-j > 0
                    temp = trackCoord(iTr).imageCoords(currentTimePoint-j,:);
                else
                    temp = NaN(1,3);
                end
                currentSpots_minusX = temp(1);
                currentSpots_minusY = temp(2);
                currentSpots_minusZ = temp(3);
                if currentTimePoint+j < timeInfo.nTimePoints
                    temp = trackCoord(iTr).imageCoords(currentTimePoint+j,:);
                else
                    temp = NaN(1,3);
                end
                currentSpots_plusX = temp(1);
                currentSpots_plusY = temp(2);
                currentSpots_plusZ = temp(3);
                
                
                trackSegTmp(toPlot-co_minus,1,1) = currentSpots_minusX;
                trackSegTmp(toPlot-co_minus,1,2) = currentSpots_minusY;
                trackSegTmp(toPlot-co_minus,1,3) = currentSpots_minusZ;
                trackSegTmp(toPlot+1,1,1) = currentSpotsX;
                trackSegTmp(toPlot+1,1,2) = currentSpotsY;
                trackSegTmp(toPlot+1,1,3) = currentSpotsZ;
                trackSegTmp(toPlot+co_minus+2,1,1) = currentSpots_plusX;
                trackSegTmp(toPlot+co_minus+2,1,2) = currentSpots_plusY;
                trackSegTmp(toPlot+co_minus+2,1,3) = currentSpots_plusZ;
                
                
            end
            
            trackSeg(:,c,:) = trackSegTmp;
            %   end
        end
        
        idx = spots(:,3) < (currentZSec+0.5) & spots(:,3) > (currentZSec-0.5);
        tracks4Plane = find(idx);
        %         spots = spots((spots(:,3) < (currentZSec+0.5) & spots(:,3) > (currentZSec-0.5)),:);
        spots = spots(idx,:);
        
        for iTr = 1:size(trackSeg,2)
            if ~(trackSeg(toPlot+1,iTr,3) < (currentZSec+1.5) && trackSeg(toPlot+1,iTr,3) > (currentZSec-1.5));
                trackSeg(:,iTr,:) = NaN;
            end
        end
        
        set(figH,'Name', sprintf('Spots %s (t = %2d / %2d, z = %2d / %2d)',imageInfo.projectName,currentTimePoint, timeInfo.nTimePoints,currentZSec,zInfo.nZSec));
        
        imshow(image,[minIn maxIn],'Border','tight');
        if ~initialPlot
            %set(gca,cameraPropertyNames,cameraProps);
            zoom reset
            set(gca,{'xlim','ylim'},L)
        end
        
        if ~isempty(spots)
            hold on
            plot(spots(:,2),spots(:,1),'g*');
            plot(trackSeg(1:toPlot+1,:,2),trackSeg(1:toPlot+1,:,1),'r')
            plot(trackSeg(toPlot+1:end,:,2),trackSeg(toPlot+1:end,:,1),'b')
            
            for iTe = 1:size(spots,1)
                textCoord = spots(iTe,1:2) + [-1 1];
                col = [0 1 0];
                sisPairNo = ceil(tracks4Plane(iTe)/2);
                if mod(tracks4Plane(iTe),2) == 0
                    sisNo = 2;
                else
                    sisNo = 1;
                end
                text(textCoord(2),textCoord(1),['\color[rgb]{',num2str(col(1)),' ',num2str(col(2)),' ',num2str(col(3)),'}','Sister Pair ID ',int2str(sis2PlotTmp(sisPairNo)),':',int2str(sisNo)])
            end
            hold off
        end
        
        % end of function plotPlanes
        
        % make visibe
        %         set(figH,'Visible','on');
    end



    function featIdx = getFeatIdx(track,noTimePoints)
        % gets the featIdx of a track for all times points, leaving a NaN
        % for no feat
        
        featIdx = NaN(noTimePoints,1);
        
        startTime = track.seqOfEvents(1,1); % get track start and end times
        endTime = track.seqOfEvents(2,1);
        
        featIdx(startTime:endTime) = track.tracksFeatIndxCG';
        
        featIdx(featIdx==0) = NaN;
        
    end


end

