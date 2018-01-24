function plotCompTrack(trackedFeatureInfo,plotX,plotY,plotA,inOneFigure)
%PLOTCOMPTRACK plots the x-coordinates, y-coordinates and/or intensities along a compound track, indicating merges, splits and gaps
%
%SYNOPSIS plotCompTrackAmp(trackedFeatureInfo,plotX,plotY,plotA,inOneFigure)
%
%INPUT  trackedFeatureInfo: Output of trackCloseGapsKalman for one track:
%                           Contains the fields:
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of 
%                              frames the compound track spans. Each row
%                              consists of 
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%       plotX             : 1 if x-coordinate is to be plotted, zero
%                           otherwise.
%       plotY             : 1 if y-coordinate is to be plotted, zero
%                           otherwise.
%       plotA             : 1 if amplitude is to be plotted, zero
%                           otherwise.
%       inOneFigure       : 1 if all plots appear in one figure window (one
%                           above the other), 0 if each figure is in its
%                           own window.
%
%OUTPUT The plot(s).
%
%REMARKS gaps are dotted black lines, splits are dash-dotted black lines
%and merges are dashed lines
%
%Khuloud Jaqaman, May 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('plotCompTrack')
    disp('--plotCompTrack: Incorrect number of input arguments!');
    return
end

%extract information from input
seqOfEvents = trackedFeatureInfo.seqOfEvents;
tracksCoordAmpCG = trackedFeatureInfo.tracksCoordAmpCG;

%get first frame, last frame and number of frames
firstFrame = seqOfEvents(1,1);
lastFrame = seqOfEvents(end,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%x-coordinates
if plotX

    %open new figure window and hold on to it
    figure, hold on
    
    %if all figures will be in the same window, specify subplot
    if inOneFigure
        subplot(plotX+plotY+plotA,1,1)
        hold on
    end

    %extract x-coordinates from input
    xCoordSequence = tracksCoordAmpCG(:,1:8:end);

    %find number of segments making compound track
    numSegments = size(xCoordSequence,1);

    %plot x-coordinates as dotted black lines, closing gaps
    for i = 1 : numSegments
        indx = find(~isnan(xCoordSequence(i,:)));
        plot(indx+firstFrame-1,xCoordSequence(i,indx),'k:');
    end

    %plot x-coordinates in color leaving gaps as blank (so that they appear as
    %dotted lines in the final figure)
    plot((firstFrame:lastFrame)',xCoordSequence','marker','.');

    %find merges and splits
    indxSplit = (find(seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4))))';
    indxMerge = (find(seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4))))';

    %go over all splits
    for iSplit = indxSplit

        %get time of splitting
        timeSplit = seqOfEvents(iSplit,1);

        %determine location in xCoordSequence
        splitLoc = timeSplit - firstFrame + 1;

        %determine index of starting track
        rowS = seqOfEvents(iSplit,3);

        %determine index of splitting track
        rowSp = seqOfEvents(iSplit,4);

        %plot split as a black dash-dotted line
        plot([timeSplit-1 timeSplit],[xCoordSequence(rowSp,splitLoc-1) ...
            xCoordSequence(rowS,splitLoc)],'k-.')

    end

    %go over all merges
    for iMerge = indxMerge

        %get time of merging
        timeMerge = seqOfEvents(iMerge,1);

        %determine location in xCoordSequence
        mergeLoc = timeMerge - firstFrame + 1;

        %determine index of ending track
        rowE = seqOfEvents(iMerge,3);

        %determine index of merging track
        rowM = seqOfEvents(iMerge,4);

        %plot merge as a black dashed line
        plot([timeMerge-1 timeMerge],[xCoordSequence(rowE,mergeLoc-1) ...
            xCoordSequence(rowM,mergeLoc)],'k--')

    end

    %put axes labels
    xlabel('frame number');
    ylabel('x-coordinate (pixels)');

    %hold off of figure if each plot is in a separate figure window or if
    %this is the last plot
    if ~inOneFigure || (~plotY && ~plotA)
        hold off
    end

end

%y-coordinates
if plotY

    %open new figure window if needed and hold on to it
    if ~inOneFigure || ~plotX
        figure, hold on
    end
    
    %if all figures will be in the same window, specify subplot
    if inOneFigure
        subplot(plotX+plotY+plotA,1,plotX+1)
        hold on
    end

    %extract y-coordinates from input
    yCoordSequence = tracksCoordAmpCG(:,2:8:end);

    %find number of segments making compound track
    numSegments = size(yCoordSequence,1);

    %plot x-coordinates as dotted black lines, closing gaps
    for i = 1 : numSegments
        indx = find(~isnan(yCoordSequence(i,:)));
        plot(indx+firstFrame-1,yCoordSequence(i,indx),'k:');
    end

    %plot x-coordinates in color leaving gaps as blank (so that they appear as
    %dotted lines in the final figure)
    plot((firstFrame:lastFrame)',yCoordSequence','marker','.');

    %find merges and splits
    indxSplit = (find(seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4))))';
    indxMerge = (find(seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4))))';

    %go over all splits
    for iSplit = indxSplit

        %get time of splitting
        timeSplit = seqOfEvents(iSplit,1);

        %determine location in yCoordSequence
        splitLoc = timeSplit - firstFrame + 1;

        %determine index of starting track
        rowS = seqOfEvents(iSplit,3);

        %determine index of splitting track
        rowSp = seqOfEvents(iSplit,4);

        %plot split as a black dash-dotted line
        plot([timeSplit-1 timeSplit],[yCoordSequence(rowSp,splitLoc-1) ...
            yCoordSequence(rowS,splitLoc)],'k-.')

    end

    %go over all merges
    for iMerge = indxMerge

        %get time of merging
        timeMerge = seqOfEvents(iMerge,1);

        %determine location in yCoordSequence
        mergeLoc = timeMerge - firstFrame + 1;

        %determine index of ending track
        rowE = seqOfEvents(iMerge,3);

        %determine index of merging track
        rowM = seqOfEvents(iMerge,4);

        %plot merge as a black dashed line
        plot([timeMerge-1 timeMerge],[yCoordSequence(rowE,mergeLoc-1) ...
            yCoordSequence(rowM,mergeLoc)],'k--')

    end

    %put axes labels
    xlabel('frame number');
    ylabel('y-coordinate (pixels)');

    %hold off of figure if each plot is in a separate figure window or if
    %this is the last plot
    if ~inOneFigure || ~plotA
        hold off
    end

end

%amplitudes
if plotA

    %open new figure window if needed and hold on to it
    if ~inOneFigure || (~plotX && ~plotY)
        figure, hold on
    end
    
    %if all figures will be in the same window, specify subplot
    if inOneFigure
        subplot(plotX+plotY+plotA,1,plotX+plotY+1)
        hold on
    end

    %extract amplitudes from input
    ampSequence = tracksCoordAmpCG(:,4:8:end);

    %find number of segments making compound track
    numSegments = size(ampSequence,1);

    %plot amplitudes as dotted black lines, closing gaps
    for i = 1 : numSegments
        indx = find(~isnan(ampSequence(i,:)));
        plot(indx+firstFrame-1,ampSequence(i,indx),'k:');
    end

    %plot amplitudes in color leaving gaps as blank (so that they appear as
    %dotted lines in the final figure)
    plot((firstFrame:lastFrame)',ampSequence','marker','.');

    %find merges and splits
    indxSplit = (find(seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4))))';
    indxMerge = (find(seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4))))';

    %go over all splits
    for iSplit = indxSplit

        %get time of splitting
        timeSplit = seqOfEvents(iSplit,1);

        %determine location in ampSequence
        splitLoc = timeSplit - firstFrame + 1;

        %determine index of starting track
        rowS = seqOfEvents(iSplit,3);

        %determine index of splitting track
        rowSp = seqOfEvents(iSplit,4);

        %plot split as a black dash-dotted line
        plot([timeSplit-1 timeSplit],[ampSequence(rowSp,splitLoc-1) ...
            ampSequence(rowS,splitLoc)],'k-.')

    end

    %go over all merges
    for iMerge = indxMerge

        %get time of merging
        timeMerge = seqOfEvents(iMerge,1);

        %determine location in ampSequence
        mergeLoc = timeMerge - firstFrame + 1;

        %determine index of ending track
        rowE = seqOfEvents(iMerge,3);

        %determine index of merging track
        rowM = seqOfEvents(iMerge,4);

        %plot merge as a black dashed line
        plot([timeMerge-1 timeMerge],[ampSequence(rowE,mergeLoc-1) ...
            ampSequence(rowM,mergeLoc)],'k--')

    end

    %put axes labels
    xlabel('frame number');
    ylabel('amplitude (normalized)');

    %hold off of figure
    hold off

end


%%%%% ~~ the end ~~ %%%%%

