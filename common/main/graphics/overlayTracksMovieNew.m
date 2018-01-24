function overlayTracksMovieNew(tracksFinal,startend,dragtailLength,...
    saveMovie,movieName,filterSigma,classifyGaps,highlightES)
%Overlays tracks obtained via trackCloseGapsKalman on movies
%
%SYNPOSIS overlayTracksMovieNew(tracksFinal,startend,dragtailLength,...
%    saveMovie,movieName,filterSigma,classifyGaps,highlightES)
%
%INPUT  tracksFinal   : Output of trackCloseGapsKalman.
%       startend      : Row vector indicating first and last frame to
%                       include in movie. Format: [startframe endframe]. 
%                       Optional. Default: [(first frame with tracks) (last frame with tracks)]            
%       dragtailLength: Length of drag tail (in frames). 
%                       Optional. Default: 10 frames.
%                       If dragtailLength = 0, then no dragtail. 
%                       To show full tracks, set dragtailLength to any 
%                       value longer than the movie.
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
%       imageRange    : Image region to make movie out of, in the form:
%                       [min pixel X, max pixel X; min pixel Y, max pixel Y].
%                       Optional. Default: Whole image.
%                       NOT IMPLEMENTED YET!!!
%
%OUTPUT If movie is to be saved, the QT movie is written into directory
%       where TIFFs are located
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
if nargin < 9;
    imageRange = [1 imSizeX; 1 imSizeY];
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

%get number of segments making each track
numSegments = zeros(numTracks,1);
for i = 1 : numTracks
    numSegments(i) = size(tracksFinal(i).tracksCoordAmpCG,1);
end

%locate the row of the first track of each compound track in the
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
    
end

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

%% make movie

%go to directory where movie will be saved
cd(dirName);

%go over all specified frames
for iFrame = 1 : size(xCoordMatAll,2)
    
    %plot image in current frame
    clf;
    imshow(ImageStack(imageRange(1,1):imageRange(1,2),...
        imageRange(2,1):imageRange(2,2),iFrame),[]); 
    hold on; 

    %show frame number
    text(10,10,num2str(iFrame),'Color','white');

    %plot tracks
    if iFrame > 1
        dragTailStart = max(iFrame-dragtailLength,1);
        xCoord2plot = (xCoordMatAll(pointStatus(:,iFrame)~=0,dragTailStart:iFrame))';
        yCoord2plot = (yCoordMatAll(pointStatus(:,iFrame)~=0,dragTailStart:iFrame))';
        plot(xCoord2plot,yCoord2plot,'Color',[1 0.7 0.7]);
    end
    
    %plot points (features + gaps)
    
    %blue stars: bad gaps
    points2plot = find(pointStatus(:,iFrame)==-2);
    plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'b*','MarkerSize',6);
    
    %cyan stars: good gaps
    points2plot = find(pointStatus(:,iFrame)==-1);
    plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'c*','MarkerSize',6);
    
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
    plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'gd','MarkerSize',6);

    %yellow diamonds: detected feature just before/after a merge
    points2plot = find(pointStatus(:,iFrame)==7);
    plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'yd','MarkerSize',6);

    %add frame to movie if movie is saved
    if saveMovie
        MakeQTMovie addfigure
    end

    %pause for a moment to see frame
    pause(0.1);

end

%finish movie
if saveMovie==1
    MakeQTMovie finish
end

%% change directory back to original
cd(startDir);

%% ~~~ end ~~~