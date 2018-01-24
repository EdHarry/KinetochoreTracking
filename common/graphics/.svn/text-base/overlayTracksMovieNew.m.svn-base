function overlayTracksMovieNew(tracksFinal,startend,dragtailLength,...
    saveMovie,movieName,filterSigma,classifyGaps,highlightES,showRaw,...
    imageRange,onlyTracks,classifyLft,diffAnalysisRes)
%OVERLAYTRACKSMOVIENEW Overlays tracks obtained via trackCloseGapsKalman on movies with variable color-coding
%
%SYNPOSIS overlayTracksMovieNew(tracksFinal,startend,dragtailLength,...
%    saveMovie,movieName,filterSigma,classifyGaps,highlightES,showRaw,...
%    imageRange,onlyTracks)
%
%INPUT  tracksFinal   : Output of trackCloseGapsKalman.
%       startend      : Row vector indicating first and last frame to
%                       include in movie. Format: [startframe endframe]. 
%                       Optional. Default: [(first frame with tracks) (last frame with tracks)]            
%       dragtailLength: Length of drag tail (in frames). 
%                       Optional. Default: 10 frames.
%                       ** If dragtailLength = 0, then no dragtail. 
%                       ** To show tracks from their beginning to their end,
%                       set dragtailLength to any value longer than the
%                       movie. 
%                       ** To show tracks statically while features dance
%                       on them, use -1.
%                       ** To show tracks from their beginning to their
%                       end, and to retain tracks even after the particle
%                       disappears, use -2.
%       saveMovie     : 1 to save movie (as Quicktime), 0 otherwise.
%                       Optional. Default: 0.
%       movieName     : filename for saving movie. 
%                       Optional. Default: TrackMovie (if saveMovie = 1).
%       filterSigma   : 0 to overlay on raw image, PSF sigma to overlay on
%                       image filtered with given filterSigma. 
%                       Optional. Default: 0.
%       classifyGaps  : 1 to classify gaps as "good" and "bad", depending
%                       on their length relative to the legnths of the
%                       segments they connect, 0 otherwise.
%                       Optional. Default: 1.
%       highlightES   : 1 to highlight track ends and starts, 0 otherwise.
%                       Optional. Default: 1.
%       showRaw       : 1 to add raw movie to the left of the movie with
%                       tracks overlaid, 2 to add raw movie at the top of
%                       the movie with tracks overlaid, 0 otherwise.
%                       Optional. Default: 0.
%       imageRange    : Image region to make movie out of, in the form:
%                       [min pixel X, max pixel X; min pixel Y, max pixel Y].
%                       Optional. Default: Whole image.
%       onlyTracks    : 1 to show only tracks without any symbols showing
%                       detections, closed gaps, merges and splits; 0 to
%                       show symbols on top of tracks.
%                       Optional. Default: 0.
%       classifyLft   : 1 to classify objects based on (1) whether they
%                       exist throughout the whole movie, (2) whether they
%                       appear OR disappear, and (3) whether they appear
%                       AND disappear; 0 otherwise. 
%                       Optional. Default: 1.
%       diffAnalysisRes: Diffusion analysis results (output of
%                       trackDiffusionAnalysis1). Needed if tracks are to be
%                       colored based on their diffusion classification.
%                       With this option, classifyGaps, highlightES and
%                       classifyLft are set to zero, regardless of input.
%                       Optional. Default: None.
%
%OUTPUT If movie is to be saved, the QT movie is written into directory
%       where TIFFs are located
%
%REMARKS Color-coding:
%        ** Without diffusion classification, all tracks have a neutral
%        color, while objects are color coded in the following way:
%               * Detected object just after appearance: Green circle.
%               * Detected object just before disappearance: Yellow
%                 circle.
%               * Detected object in middle of trajectory that spans
%                 whole movie: White circle.
%               * Detected object in middle of trajectory that appears OR
%                 disappears within movie: Magenta circle.
%               * Detected object in middle of trajectory that appears AND
%                 disappears within movie: Red circle.
%               * Gap that is short than both segments it connects: Cyan
%                 star.
%               * Gap that is longer than at least one ofthe segments it
%                 connects: Blue star.
%               * Object before and after splitting: Green diamond.
%               * OBject before and after merging: Yellow diamond.
%           When classifyGaps = 0, all gaps are cyan.
%           When highlighES = 0, no green and yellow circles.
%           When classifyLft = 0, all objets in middle of trajectory are white.
%
%       ** With diffusion classification, all objects and gaps have neutral
%       color (merges and splits are diamonds), while tracks are
%       color coded in the following way:
%               * Type 1: Linear + 1D confined diffusion: Orange.
%               * Type 2: Linear + 1D normal diffusion: Red.
%               * Type 3: Linear + 1D super diffusion: Green.
%               * Type 4: Linear + too short for 1D classification: Yellow.
%               * Type 5: Random/Unclassified + 2D confined diffusion: Blue.
%               * Type 6: Random/Unclassified + 2D normal diffusion: Cyan.
%               * Type 7: Random/Unclassified + 2D super diffusion: Magenta.
%               * Type 8: Random + too short for 2D classification: Purple.
%               * Type 0: Too short for any analysis: Light pink.
%
%Khuloud Jaqaman, August 2007

%% input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--overlayTracksMovieNew: Incorrect number of input arguments!');
    return
end

%get first and last frames where there are tracks
allEvents = vertcat(tracksFinal.seqOfEvents);
tracksFirstFrame = min(allEvents(:,1));
tracksLastFrame = max(allEvents(:,1));

%record directory before start of function
startDir = pwd;

%ask user for images
[fName,dirName] = uigetfile('*.tif','specify first image in the stack - specify very first image, even if not to be plotted');

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileList = getFileStackNames([dirName,fName]);
    numFrames = length(outFileList);
    
    %read first image to get image size
    currentImage = imread(outFileList{1});
    [isx,isy] = size(currentImage);

else %else, exit
    
    disp('--overlayTracksMovieNew: Bad file selection');
    return

end

%check startend and assign default if necessary
if nargin < 2 || isempty(startend)
    startend = [tracksFirstFrame tracksLastFrame];
else
    tracksFirstFrame = min(tracksFirstFrame,startend(1));
    tracksLastFrame = max(tracksLastFrame,startend(2));
end

%check dragtailLength and assign default if not necessary
if nargin < 3 || isempty(dragtailLength)
    dragtailLength = 10;
end

%check whether to save movie
if nargin < 4 || isempty(saveMovie)
    saveMovie = 0;
end

%check name for saving movie
if saveMovie && (nargin < 5 || isempty(movieName))
    movieName = 'trackMovie.mov';
end

%check whether to use filtered images
if nargin < 6 || isempty(filterSigma)
    filterSigma = 0;
end

%check whether to color-code gaps
if nargin < 7 || isempty(classifyGaps)
    classifyGaps = 1;
end

%check whether to highligh track starts and ends
if nargin < 8 || isempty(highlightES)
    highlightES = 1;
end

%check whether to put raw movie adjacent to movie with tracks overlaid
if nargin < 9 || isempty(showRaw)
    showRaw = 0;
end

%keep only the frames of interest
outFileList = outFileList(startend(1):startend(2));

%read specified images into ImageStack
ImageStack = zeros(isx,isy,length(outFileList));
h = waitbar(0,'reading images');
for i=1:length(outFileList)
    currentImage = imread(outFileList{i});
    ImageStack(:,:,i) = currentImage;
    if (mod(i,10)==0), waitbar(i/length(outFileList)); end
end
close(h);

%get image size
[imSizeX,imSizeY] = size(ImageStack(:,:,1));

%check whether an area of interest was input
if nargin < 10 || isempty(imageRange)
    imageRange = [1 imSizeX; 1 imSizeY];
end

%check whether to plot tracks only or also symbols
if nargin < 11 || isempty(onlyTracks)
    onlyTracks = 0;
end

%check whether to classify lifetime
if nargin < 12 || isempty(classifyLft)
    classifyLft = 1;
end

%check whether to color-code tracks based on diffusion classification
if nargin < 13 || isempty(diffAnalysisRes)
    diffAnalysisRes = [];
else
    classifyGaps = 0;
    highlightES = 0;
    classifyLft = 0;
end

%filter images if requested
if filterSigma
    for i = 1 : length(outFileList)
        imageStack(:,:,i) = Gauss2D(imageStack(:,:,i),filterSigma);
    end
end

%initialize QT movie if it is to be saved
if saveMovie
    evalString = ['MakeQTMovie start ',movieName];
    eval(evalString);
end

%% store track positions, get track status and point status

%get number of tracks
numTracks = length(tracksFinal);

%get track start and end times
trackSEL = getTrackSEL(tracksFinal);

%give tracks status based on the frames they span: 
%2: track exists throughout movie
%1: track exists either in first frame or in last frame
%0: track does not exist in both first frame and last frame
trackStatus  = (trackSEL(:,1) == tracksFirstFrame) + (trackSEL(:,2) == tracksLastFrame);

%give all tracks same classification if lifetime classification not
%requested
if classifyLft == 0
    trackStatus(:) = 2;
end

%get number of segments making each track
numSegments = zeros(numTracks,1);
for i = 1 : numTracks
    numSegments(i) = size(tracksFinal(i).tracksCoordAmpCG,1);
end

%locate the row of the first segment of each compound track in the
%big matrices of all tracks (to be constructed in the next step)
trackStartRow = ones(numTracks,1);
for iTrack = 2 : numTracks
    trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
end

%find total number of segments in all tracks (i.e. number of rows in big
%matrices)
numSegmentsTracks = trackStartRow(end)+numSegments(end)-1;

%construct a matrix indicating point status in big matrices:
%-2: bad gap (gap length > either segment length on its sides)
%-1: good gap (gap length < both segment lengths on its sides)
%0 : before track start or after track end
%1 : detected feature in the middle of a track of trackStatus = 0
%2 : detected feature in the middle of a track of trackStatus = 1
%3 : detected feature in the middle of a track of trackStatus = 2
%4 : detected feature just after a birth
%5 : detected feature just before a death
%6 : detected feature just after a split
%7 : detected feature just before a merge
pointStatus = zeros(numSegmentsTracks,numFrames);

%put all tracks together in one big matrix
%put the x-coordinates in one matrix and the y-coordinates in another
%indicate the status of each point
xCoordMatAll = NaN*ones(numSegmentsTracks,numFrames);
yCoordMatAll = xCoordMatAll;
for iTrack = 1 : numTracks
    
    %get track start and end times
    startTime = trackSEL(iTrack,1);
    endTime   = trackSEL(iTrack,2);
    
    %store x-coordinates
    xCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = ...
        tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end);
    
    %store y-coordinates
    yCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = ...
        tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end);
    
    %assign point status for features in the middle of the track
    pointStatus(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = trackStatus(iTrack) + 1;
    
    %get sequence of events of track
    seqOfEvents = tracksFinal(iTrack).seqOfEvents;
    
    if highlightES

        %assign point status for features just after a birth
        points2consider = find(seqOfEvents(:,2)==1 & isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1)~=tracksFirstFrame)';
        for iPoint = points2consider
            pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
                seqOfEvents(iPoint,1)) = 4;
        end

        %assign point status for features just before a death
        points2consider = find(seqOfEvents(:,2)==2 & isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1)~=tracksLastFrame)';
        for iPoint = points2consider
            pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
                seqOfEvents(iPoint,1)) = 5;
        end
        
    end
    
    %assign point status for features just after and before a split
    %also, in the frame just before splitting, give the
    %splitting track the position of the track it split from 
    points2consider = find(seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)))';
    for iPoint = points2consider
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 6;
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,4)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 6;
        xCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1) = xCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1);
        yCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1) = yCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1);
    end
    
    %assign point status for features just before and after a merge
    %also, in the frame just after merging, give the
    %merging track the position of the track it merged from 
    points2consider = find(seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)))';
    for iPoint = points2consider
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 7;
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,4)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 7;
        xCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)) = xCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1));
        yCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)) = yCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1));
    end
    
end %(for iTrack = 1 : numTracks)

%find gaps in tracks
gapInfo = findTrackGaps(tracksFinal);

%for gaps, assign the position as that of the feature before the gap
%also, assign a point status of -1
for iGap = 1 : size(gapInfo,1)
    
    iTrack = gapInfo(iGap,1);
    iSegment = gapInfo(iGap,2);
    iStart = gapInfo(iGap,3);
    gapLength = gapInfo(iGap,4);
    if classifyGaps
        if gapInfo(iGap,5) <= 1 && gapInfo(iGap,6) <= 1
            gapType = -1;
        else
            gapType = -2;
        end
    else
        gapType = -1;
    end
    
    xCoordMatAll(trackStartRow(iTrack)+iSegment-1,iStart:iStart+gapLength-1) = ...
        xCoordMatAll(trackStartRow(iTrack)+iSegment-1,iStart-1); %x-coordinates
    
    yCoordMatAll(trackStartRow(iTrack)+iSegment-1,iStart:iStart+gapLength-1) = ...
        yCoordMatAll(trackStartRow(iTrack)+iSegment-1,iStart-1); %y-coordinates
    
    pointStatus(trackStartRow(iTrack)+iSegment-1,iStart:iStart+gapLength-1) = gapType; %point status
    
end

%retain in the big matrices only the frames of interest
xCoordMatAll = xCoordMatAll(:,startend(1):startend(2));
yCoordMatAll = yCoordMatAll(:,startend(1):startend(2));

%% get track classification based on diffusion analysis

%assign default classification as 0, which means undetermined, and which is
%the only classification when there are no diffusion analysis results
trackClass = zeros(numSegmentsTracks,1);

%if there are diffusion analysis results ...
if ~isempty(diffAnalysisRes)

    %get track segment types from diffusion analysis
    trackSegmentType = vertcat(diffAnalysisRes.classification);

    %assign classes
    trackClass(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 1) = 1; %linear + 1D confined (orange)
    trackClass(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 2) = 2; %linear + 1D normal (red)
    trackClass(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 3) = 3; %linear + 1D super (green)
    trackClass(trackSegmentType(:,1) == 1 & isnan(trackSegmentType(:,3))) = 4; %linear + too short (yellow)
    trackClass(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 1) = 5; %random/unclassified + 2D confined (blue)
    trackClass(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 2) = 6; %random/unclassified + 2D normal (cyan)
    trackClass(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 3) = 7; %random/unclassified + 2D super (magenta)
    trackClass(trackSegmentType(:,1) == 0 & isnan(trackSegmentType(:,2))) = 8; %random + too short (purple)    
    
end

%% make movie

%go to directory where movie will be saved
cd(dirName);

%go over all specified frames
for iFrame = 1 : size(xCoordMatAll,2)
    
    %plot image in current frame and show frame number
    clf;
    switch showRaw
        case 1
            axes('Position',[0 0 0.495 1]);
            imshow(ImageStack(:,:,iFrame),[]);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
                textDeltaCoord,num2str(iFrame),'Color','white');
            axes('Position',[0.505 0 0.495 1]);
            imshow(ImageStack(:,:,iFrame),[]);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
        case 2
            axes('Position',[0 0.505 1 0.495]);
            imshow(ImageStack(:,:,iFrame),[]);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
                textDeltaCoord,num2str(iFrame),'Color','white');
            axes('Position',[0 0 1 0.495]);
            imshow(ImageStack(:,:,iFrame),[]);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
        otherwise
            axes('Position',[0 0 1 1]);
            imshow(ImageStack(:,:,iFrame),[]);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
                textDeltaCoord,num2str(iFrame),'Color','white');
    end

    %get tracks to plot
    plotOrNot = 0;
    if dragtailLength >= 0 %plot tracks dynamically
        
        if iFrame > 1
            dragTailStart = max(iFrame-dragtailLength,1);
            indx2keep = find(pointStatus(:,iFrame)~=0);
            xCoord2plot = (xCoordMatAll(indx2keep,dragTailStart:iFrame))';
            yCoord2plot = (yCoordMatAll(indx2keep,dragTailStart:iFrame))';
            trackClass2plot = trackClass(indx2keep);
            plotOrNot = 1;
        end
        
    elseif dragtailLength == -1 %plot tracks statically
        
        xCoord2plot = (xCoordMatAll)';
        yCoord2plot = (yCoordMatAll)';
        trackClass2plot = trackClass;
        plotOrNot = 1;
        
    elseif dragtailLength == -2 %plot tracks dynamically but keep them after they disappear
        
        if iFrame > 1
            xCoord2plot = xCoordMatAll(:,1:iFrame)';
            yCoord2plot = yCoordMatAll(:,1:iFrame)';
            trackClass2plot = trackClass;
            plotOrNot = 1;
        end
        
    end

    %plot tracks
    %color-code dragtail based on diffusion analysis if supplied
    if plotOrNot
        plot(xCoord2plot(:,trackClass2plot==0),yCoord2plot(:,trackClass2plot==0),...
            'Color',[1 0.7 0.7],'LineWidth',2); %light pink
        plot(xCoord2plot(:,trackClass2plot==1),yCoord2plot(:,trackClass2plot==1),...
            'Color',[1 0.7 0],'LineWidth',2); %orange
        plot(xCoord2plot(:,trackClass2plot==2),yCoord2plot(:,trackClass2plot==2),...
            'Color','r','LineWidth',2); %[1 0 0]
        plot(xCoord2plot(:,trackClass2plot==3),yCoord2plot(:,trackClass2plot==3),...
            'Color','g','LineWidth',2); %[0 1 0]
        plot(xCoord2plot(:,trackClass2plot==4),yCoord2plot(:,trackClass2plot==4),...
            'Color','y','LineWidth',2); %[1 1 0]
        plot(xCoord2plot(:,trackClass2plot==5),yCoord2plot(:,trackClass2plot==5),...
            'Color','b','LineWidth',2); %[0 0 1]
        plot(xCoord2plot(:,trackClass2plot==6),yCoord2plot(:,trackClass2plot==6),...
            'Color','c','LineWidth',2); %[0 1 1]
        plot(xCoord2plot(:,trackClass2plot==7),yCoord2plot(:,trackClass2plot==7),...
            'Color','m','LineWidth',2); %[1 0 1]
        plot(xCoord2plot(:,trackClass2plot==8),yCoord2plot(:,trackClass2plot==8),...
            'Color',[0.6 0 1],'LineWidth',2); %purple
    end

    %plot points (features + gaps + merges + splits)
    if ~onlyTracks

        %blue stars: bad gaps
        points2plot = find(pointStatus(:,iFrame)==-2);
        plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'b*','MarkerSize',6);

        %cyan stars: good gaps
        points2plot = find(pointStatus(:,iFrame)==-1);
        if isempty(diffAnalysisRes)
            plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'c*','MarkerSize',6);
        else
            plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'w*','MarkerSize',6);
        end

        %red circles: detected feature in the middle of track with status 0
        points2plot = find(pointStatus(:,iFrame)==1);
        plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'ro','MarkerSize',5);

        %magenta circles: detected feature in the middle of track with status 1
        points2plot = find(pointStatus(:,iFrame)==2);
        plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'mo','MarkerSize',5);

        %white circles: detected feature in the middle of track with status 2
        points2plot = find(pointStatus(:,iFrame)==3);
        plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'wo','MarkerSize',5);

        %green circles: detected feature just after birth
        points2plot = find(pointStatus(:,iFrame)==4);
        plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'go','MarkerSize',5);

        %yellow circles: detected feature just before death
        points2plot = find(pointStatus(:,iFrame)==5);
        plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'yo','MarkerSize',5);

        %green diamonds: detected feature just before/after a split
        points2plot = find(pointStatus(:,iFrame)==6);
        if isempty(diffAnalysisRes)
            plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'gd','MarkerSize',6);
        else
            plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'wd','MarkerSize',6);
        end

        %yellow diamonds: detected feature just before/after a merge
        points2plot = find(pointStatus(:,iFrame)==7);
        if isempty(diffAnalysisRes)
            plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'yd','MarkerSize',6);
        else
            plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'wd','MarkerSize',6);
        end
            
    end

    %add frame to movie if movie is saved
    if saveMovie
        MakeQTMovie addfigure
    end

    %pause for a moment to see frame
    pause(0.1);

end %(for iFrame = 1 : size(xCoordMatAll,2))

%finish movie
if saveMovie==1
    MakeQTMovie finish
end

%% change directory back to original
cd(startDir);

%% ~~~ end ~~~