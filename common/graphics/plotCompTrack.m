function plotCompTrack(trackedFeatureInfo,plotX,plotY,plotA,inOneFigure,...
    plotAggregState)
%PLOTCOMPTRACK plots the x-coordinates, y-coordinates and/or intensities along a compound track, indicating merges, splits and gaps
%
%SYNOPSIS plotCompTrackAmp(trackedFeatureInfo,plotX,plotY,plotA,inOneFigure,...
%    plotAggregState)
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
%           .aggregationState: This field results from running the function
%                              aggregStateFromCompTracks. Only needed if
%                              aggregation state is to be plotted.
%       plotX             : 1 if x-coordinate is to be plotted, zero
%                           otherwise. Optional. Default: 1.
%       plotY             : 1 if y-coordinate is to be plotted, zero
%                           otherwise. Optional. Default: 1.
%       plotA             : 1 if amplitude is to be plotted, zero
%                           otherwise. Optional. Default: 1.
%       inOneFigure       : 1 if all plots appear in one figure window (one
%                           above the other), 0 if each figure is in its
%                           own window. Optional. Default: 1.
%       plotAggregState   : 1 to plot particle aggregation state (if
%                           supplied), 0 otherwise. Optional. Default: 1.
%
%OUTPUT The plot(s).
%
%REMARKS gaps are dotted black lines, splits are dash-dotted black lines
%and merges are dashed lines
%
%Khuloud Jaqaman, May 2007

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--plotCompTrack: Incorrect number of input arguments!');
    return
end

%assign defaults if parameters were not input
if nargin < 2 || isempty(plotX)
    plotX = 1;
end
if nargin < 3 || isempty(plotY)
    plotY = 1;
end
if nargin < 4 || isempty(plotA)
    plotA = 1;
end
if nargin < 5 || isempty(inOneFigure)
    inOneFigure = 1;
end
if nargin < 6 || isempty(plotAggregState)
    plotAggregState = 1;
end

%extract information from input
seqOfEvents = trackedFeatureInfo.seqOfEvents;
tracksCoordAmpCG = trackedFeatureInfo.tracksCoordAmpCG;
if isfield(trackedFeatureInfo,'aggregState')
    aggregState = trackedFeatureInfo.aggregState;
else
    plotAggregState = 0;
end

%determine whether track is sampled regularly or with doubled frequency
doubleFreq = mod(seqOfEvents(1,1)*2,2)==1;

%% Plotting

%x-coordinates
if plotX

    %open new figure window and hold on to it
    figure, hold on

    %if all figures will be in the same window, specify subplot
    if inOneFigure
        subplot(plotX+plotY+plotA+plotAggregState,1,1)
        hold on
    end

    %extract x-coordinates from input
    xCoordSequence = tracksCoordAmpCG(:,1:8:end);
    
    %plot the x-coordinate, taking into account closed gaps, merges and
    %splits and sampling frequency doubling
    plotCompTrackCore(xCoordSequence,seqOfEvents,doubleFreq);

    %put axes labels
    xlabel('Frame number');
    ylabel('X-coordinate (pixels)');

    %hold off of figure if each plot is in a separate figure window or if
    %this is the last plot
    if ~inOneFigure || (~plotY && ~plotA && ~plotAggregState)
        hold off
    end

end %(if plotX)

%y-coordinates
if plotY

    %open new figure window if needed and hold on to it
    if ~inOneFigure || ~plotX
        figure, hold on
    end

    %if all figures will be in the same window, specify subplot
    if inOneFigure
        subplot(plotX+plotY+plotA+plotAggregState,1,plotX+1)
        hold on
    end

    %extract y-coordinates from input
    yCoordSequence = tracksCoordAmpCG(:,2:8:end);

    %plot the y-coordinate, taking into account closed gaps, merges and
    %splits and sampling frequency doubling
    plotCompTrackCore(yCoordSequence,seqOfEvents,doubleFreq);

    %put axes labels
    xlabel('Frame number');
    ylabel('Y-coordinate (pixels)');

    %hold off of figure if each plot is in a separate figure window or if
    %this is the last plot
    if ~inOneFigure || (~plotA && ~plotAggregState)
        hold off
    end

end %(if plotY)

%amplitudes
if plotA

    %open new figure window if needed and hold on to it
    if ~inOneFigure || (~plotX && ~plotY)
        figure, hold on
    end

    %if all figures will be in the same window, specify subplot
    if inOneFigure
        subplot(plotX+plotY+plotA+plotAggregState,1,plotX+plotY+1)
        hold on
    end

    %extract amplitudes from input
    ampSequence = tracksCoordAmpCG(:,4:8:end);

    %plot the amplitudes, taking into account closed gaps, merges and
    %splits and sampling frequency doubling
    plotCompTrackCore(ampSequence,seqOfEvents,doubleFreq);

    %put axes labels
    xlabel('Frame number');
    ylabel('Amplitude (a.u.)');

    %hold off of figure if each plot is in a separate figure window or if
    %this is the last plot
    if ~inOneFigure || (~plotAggregState)
        hold off
    end

end %(if plotA)

%aggregation state
if plotAggregState

    %open new figure window if needed and hold on to it
    if ~inOneFigure || (~plotX && ~plotY && ~plotA)
        figure, hold on
    end

    %if all figures will be in the same window, specify subplot
    if inOneFigure
        subplot(plotX+plotY+plotA+plotAggregState,1,plotX+plotY+plotA+1)
        hold on
    end

    %replace zeros with NaNs in aggregation state matrix
    aggregState(aggregState==0) = NaN;
    
    %plot the aggregation states, taking into account closed gaps, merges and
    %splits and sampling frequency doubling
    plotCompTrackCore(aggregState,seqOfEvents,doubleFreq);

    %put axes labels
    xlabel('Frame number');
    ylabel('Aggregation state');

    %hold off of figure
    hold off

end %(if plotAggregState)

%% Subfunction

function plotCompTrackCore(valuesMatrix,seqOfEvents,doubleFreq)

%get first frame, last frame and number of frames
firstFrame = seqOfEvents(1,1);
lastFrame = seqOfEvents(end,1);

%get sampling frequency
samplingFreq = 1 / (1+doubleFreq);

%get number of segments making compound track
numSegments = size(valuesMatrix,1);

%plot values as dotted black lines, closing gaps
for i = 1 : numSegments
    indx = find(~isnan(valuesMatrix(i,:)));
    plot((indx-1)*samplingFreq+firstFrame,valuesMatrix(i,indx),'k:');
end

%plot values in color, leaving gaps as blank (so that they appear as
%dotted lines in the final figure)
plot((firstFrame:samplingFreq:lastFrame)',valuesMatrix','marker','.');

%find merges and splits
indxSplit = (find(seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4))))';
indxMerge = (find(seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4))))';

%go over all splits
for iSplit = indxSplit

    %get time of splitting
    timeSplit = seqOfEvents(iSplit,1);

    %determine location in valuesMatrix
    splitLoc = (timeSplit - firstFrame) / samplingFreq + 1;

    %determine index of starting track
    rowS = seqOfEvents(iSplit,3);

    %determine index of splitting track
    rowSp = seqOfEvents(iSplit,4);

    %plot split as a black dash-dotted line
    plot([timeSplit-samplingFreq timeSplit],[valuesMatrix(rowSp,splitLoc-1) ...
        valuesMatrix(rowS,splitLoc)],'k-.')

end

%go over all merges
for iMerge = indxMerge

    %get time of merging
    timeMerge = seqOfEvents(iMerge,1);

    %determine location in valuesMatrix
    mergeLoc = (timeMerge - firstFrame) / samplingFreq + 1;

    %determine index of ending track
    rowE = seqOfEvents(iMerge,3);

    %determine index of merging track
    rowM = seqOfEvents(iMerge,4);

    %plot merge as a black dashed line
    plot([timeMerge-samplingFreq timeMerge],[valuesMatrix(rowE,mergeLoc-1) ...
        valuesMatrix(rowM,mergeLoc)],'k--')

end


%% ~~~ the end ~~~

