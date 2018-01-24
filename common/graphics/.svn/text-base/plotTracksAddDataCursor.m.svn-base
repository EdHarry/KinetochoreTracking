function dcmObj = plotTracksAddDataCursor(figureHandle,tracksFinal,plotOpt,flipXY,moreData)
%PLTOTRACKSADDDATACURSOR adds a data cursor to a plot of tracks
%
% SYNOPSIS: dcmObj = plotTracksAddDataCursor(figureHandle,tracksFinal,plotOpt,flipXY,moreData)
%
% INPUT figureHandle : handle to figure to which dataCursor should be added
%		tracksFinal: tracksFinal structure (output of tracksCloseGapsKalman.m)
%		plotOpt: (opt) what to plot. Default: 1
%           1: basic list (+ moreData if nonempty)
%     --TBD 2: basic list + distances to other points (+ moreData if nonempty)
%		flipXY: (opt) true if x and y are flipped. Default: false
%		moreData: (opt) structure of length nTimepoints with any fields of
%                 length nFeatures(t). The content of the structure will be
%                 plotted along with the data.
%
% OUTPUT dcmObj : datacursormanager-object
%
% REMARKS
%
% created with MATLAB ver.: 7.8.0.8205 (R2009a) Beta (Mac Intel 64-bit) on Mac OS X  Version: 10.5.6 Build: 9G2141
%
% created by: jonas
% DATE: 14-Apr-2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check input
if nargin < 2 || isempty(tracksFinal) || isempty(figureHandle)
    error('plotTracksAddDataCursor needs a nonempty figure handle and a nonempty tracksFinal structure')
end

if nargin < 3 || isempty(plotOpt)
    plotOpt = 1;
elseif ~isscalar(plotOpt)
    error('plotOpt needs to be a numeric scalar')
end

if nargin < 4 || isempty(flipXY)
    flipXY = false;
elseif ~isscalar(flipXY)
    error('flipXY needs to be a numeric scalar')
end

if nargin < 5 || isempty(moreData)
    moreData = [];
elseif ~isstruct(moreData)
    error('moreData needs to be a structure of length nTimepoints with fields that contain an entry for each feature')
end

%% Add data Cursor to figure
dcmObj = datacursormode(figureHandle);

% add the callback
dcmObj.updateFcn = @(u,v)(plotTracksDataCursorUpdateFcn(u,v,tracksFinal,plotOpt,flipXY,moreData,dcmObj));

end % main function


%% plotTracksDataCursorUpdateFcn
function txt = plotTracksDataCursorUpdateFcn(empt,evt_obj,tracksFinal,plotOpt,flipXY,moreData,dcmObj)
%PLOTTRACKSDATACURSORUPDATEFCN is the update function for the data cursor object for plotting tracks
% see main help for input descriptions.
% empt, evt_obj are reserved by Matlab


%% Read minimum info

% get coords
coords = zeros(1,3);
pos = get(evt_obj,'Position');
coords(1:length(pos)) = pos;
if flipXY
    coords = coords([2,1,3]);
end

% find coords in tracksFinal by checking whether coords exist inside any
% track. This may fail for pixel coords and gridded data, because one
% position may show up multiple times. Fix it once it is necessary
nTracks = length(tracksFinal);
for iTrack = 1:nTracks
    [segment,timeIdx] = find(tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end)==coords(1) & ...
        tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end)==coords(2) & ...
        tracksFinal(iTrack).tracksCoordAmpCG(:,3:8:end)==coords(3));
    if ~isempty(segment)
        % shift time so that it agrees with timepoints in moreData
        timePoint = timeIdx + tracksFinal(iTrack).seqOfEvents(1) - 1;
        break % stop looking any further
    end
end

if isempty(segment)
    error('no correspoinding track found in data. Check for flipped XY')
end

% Get amp
amp = tracksFinal(iTrack).tracksCoordAmpCG(segment,(timeIdx-1)*8+4);

% find segment start/end
segStart = tracksFinal(iTrack).seqOfEvents(...
    tracksFinal(iTrack).seqOfEvents(:,2) == 1 & ...
    tracksFinal(iTrack).seqOfEvents(:,3) == segment,1);
segEnd = tracksFinal(iTrack).seqOfEvents(...
    tracksFinal(iTrack).seqOfEvents(:,2) == 2 & ...
    tracksFinal(iTrack).seqOfEvents(:,3) == segment,1);

featIdx = tracksFinal(iTrack).tracksFeatIndxCG(segment,timeIdx);

% make basic text
switch plotOpt
    case 1
        txt = [basicText;moreText];
    case 2
        % add distance text later
        distText = {};
        txt = [basicText;distText;moreText];
    otherwise
        error('plotOpt not implemented yet')
end

%% nested functions
    function txt = basicText
        txt = {sprintf('X/Y/Z/T : %6.2f/%6.2f/%6.2f/%i',coords,timePoint);...
            sprintf('Tr/Seg/SSt/SE/idx : %i/%i/%i/%i/%i',iTrack,segment,segStart,segEnd,featIdx);...
            sprintf('amp : %6.2f',amp);...
            };
    end

    function txt = moreText
        % check whether there is moreData
        txt = {};
        if ~isempty(moreData)
            fn = fieldnames(moreData);
            for name = fn(:)'
                name = name{1};
                txt{end+1,1} = ...
                    sprintf('%s : %6.2f',name,moreData(timePoint).(name)(featIdx));
            end
        end
    end
end
