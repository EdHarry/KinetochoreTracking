function dataStruct = makiFindAllPoles_oldVersion( dataStruct )
%MAKIFINDALLPOLES Will continually update poles tracks until all are found
%(using the new pole_tracks)
% EHarry Jan 2012

dataStruct = makiFindPoles(dataStruct); % first find the poles (1st guess)
dataStruct = makiPoleReferenceFrame(dataStruct); % then trasform the coords
dataStruct = makiUpdateClassPole(dataStruct); % update phase and class
dataStruct = makiGenerateTracks_pole(dataStruct); % then track

[newPoleTracks,found] = makiCheckForMorePoleTracks( dataStruct ); % check for new poles

while found
    dataStruct.poles = newPoleTracks; % update the poles
    dataStruct.poleReferenceFrame = [];
    dataStruct = makiPoleReferenceFrame(dataStruct); % re-transform
    dataStruct.planeFit_pole = [];
    dataStruct = makiUpdateClassPole(dataStruct); % re-update
    dataStruct.tracks_pole = [];
    dataStruct = makiGenerateTracks_pole(dataStruct); % re-track
    [newPoleTracks,found] = makiCheckForMorePoleTracks( dataStruct ); % re-check
end

% edit, EHarry April 2012, we could potentially now find poles that are
% actually kinetochores or vica versa, check for this and use this
% stragegy if a conflict is found:
%
% if a pole has its centriol then remove the false kinetochores (this is because the two centriols will have a high xcorr and low av. dis
% meaning that they are most probably poles and not kinetochores)
%
% if a pole does not have its centriol then remove the poles, this is
% because if were really a pole then the other centriol (the other sister kinetochore) would have been
% found in the centriol search
%
% remeber that if the sisterList is edited then redo the update class and
% fameAlignment


% end if no poles
if isempty(dataStruct.poles)
    return
end

% get idxs of pole tracks, let's not worry about poles from the pole-based
% tracks, if there are some then there's not going to be a sister
% kinetochore pair there
pole1TrackIdx = dataStruct.poles.pole1TrackIdx;
pole2TrackIdx = dataStruct.poles.pole2TrackIdx;

% get sister tracks idxs
sisTrackIdx = dataStruct.sisterList(1).trackPairs(:,1:2);

% check for overlaps
pole1SisOverlap = ismember(pole1TrackIdx,sisTrackIdx(:));
pole2SisOverlap = ismember(pole2TrackIdx,sisTrackIdx(:));

% if no overlaps then end
if ~any([pole1SisOverlap pole2SisOverlap])
    return
end

% if two poles (centriols) then remove the false kinetochores
falseSis1 = [];
falseSis2 = [];
if length(pole1SisOverlap) == 2
    falseSis1 = [find(pole1TrackIdx(1)==sisTrackIdx(:,1)) find(pole1TrackIdx(1)==sisTrackIdx(:,2)) find(pole1TrackIdx(2)==sisTrackIdx(:,1)) find(pole1TrackIdx(2)==sisTrackIdx(:,2))]; % check both columns
    falseSis1 = unique(falseSis1);
end
if length(pole2SisOverlap) == 2
    falseSis2 = [find(pole2TrackIdx(1)==sisTrackIdx(:,1)) find(pole2TrackIdx(1)==sisTrackIdx(:,2)) find(pole2TrackIdx(2)==sisTrackIdx(:,1)) find(pole2TrackIdx(2)==sisTrackIdx(:,2))]; % check both columns
    falseSis2 = unique(falseSis2);
end

if ~isempty(falseSis1) || ~isempty(falseSis2)
    falseSis = [falseSis1 falseSis2];
else
    falseSis = [];
end

if ~isempty(falseSis)
    % copy out the sisterList and trackPairs
    sisterList = dataStruct.sisterList;
    trackPairs = dataStruct.sisterList(1).trackPairs;
    
    % edit out the false sister(s)
    sisterList(falseSis) = [];
    trackPairs(falseSis,:) = [];
    
    % put back
    dataStruct.sisterList = sisterList;
    dataStruct.sisterList(1).trackPairs = trackPairs;
    
    % redo updateClass and frameAlignment
    dataStruct = makiUpdateClass(dataStruct);
    dataStruct = makiAlignFrames(dataStruct);
end



% if only one pole then remove the poles
if length(pole1SisOverlap) == 1
    if pole1SisOverlap
        % delete all the false pole info
        dataStruct.poles = [];
        dataStruct.poleReferenceFrame = [];
        dataStruct.planeFit_pole = [];
        dataStruct.tracks_pole = [];
        % we can end here now
        return
    end
end
if length(pole2SisOverlap) == 1
    if pole2SisOverlap
        % delete all the false pole info
        dataStruct.poles = [];
        dataStruct.poleReferenceFrame = [];
        dataStruct.planeFit_pole = [];
        dataStruct.tracks_pole = [];
        % we can end here now
        return
    end
end

end

