function alignedPoleCoords = makiGetAlignedPoleCoords( dataStruct )
%MAKIGETALIGNEDPOLECOORDS Gets the aligned (plate based) coordinates of the
%spindle poles from a makiPole project
%   EHarry Feb 2012


% get pole struct, quit if there isn't one
if ~isfield(dataStruct,'poles') || isempty(dataStruct.poles)
    alignedPoleCoords = [];
    return
end
poles = dataStruct.poles;

% get the plate based tracks
if ~isfield(dataStruct,'tracks') || isempty(dataStruct.tracks)
    alignedPoleCoords = [];
    return
end
tracks = dataStruct.tracks;

% get the pole track indexes
pole1Idx = poles.pole1TrackIdx;
pole2Idx = poles.pole2TrackIdx;

% if any of the indexes are -ve then the tracks come from the pole based
% tracks
if any(pole1Idx < 0)
    getPoleTracks1 = 1;
else
    getPoleTracks1 = 0;
end
if any(pole2Idx < 0)
    getPoleTracks2 = 1;
else
    getPoleTracks2 = 0;
end

% get the pole tracks if needed
if any([getPoleTracks1 getPoleTracks2])
    if ~isfield(dataStruct,'tracks_pole') || isempty(dataStruct.tracks_pole)
        alignedPoleCoords = [];
        return
    end
    tracks_pole = dataStruct.tracks_pole;
end


% get the featIdx, start and end of the tracks
featIdxTemp = zeros(1,dataStruct.dataProperties.movieSize(4));
c = 0;
for iPole1 = pole1Idx
    if iPole1 < 0
        c = c + 1;
        pole1Feats(c,:) = featIdxTemp;
        pole1Feats(c,tracks_pole(-iPole1).seqOfEvents(1,1):tracks_pole(-iPole1).seqOfEvents(2,1)) = tracks_pole(-iPole1).tracksFeatIndxCG;
    else
        c = c + 1;
        pole1Feats(c,:) = featIdxTemp;
        pole1Feats(c,tracks(iPole1).seqOfEvents(1,1):tracks(iPole1).seqOfEvents(2,1)) = tracks(iPole1).tracksFeatIndxCG;
    end
end
c = 0;
for iPole2 = pole2Idx
    if iPole2 < 0
        c = c + 1;
        pole2Feats(c,:) = featIdxTemp;
        pole2Feats(c,tracks_pole(-iPole2).seqOfEvents(1,1):tracks_pole(-iPole2).seqOfEvents(2,1)) = tracks_pole(-iPole2).tracksFeatIndxCG;
    else
        c = c + 1;
        pole2Feats(c,:) = featIdxTemp;
        pole2Feats(c,tracks(iPole2).seqOfEvents(1,1):tracks(iPole2).seqOfEvents(2,1)) = tracks(iPole2).tracksFeatIndxCG;
    end
end


% get the aligned coordinates from the frameAlignment struct
if ~isfield(dataStruct,'frameAlignment') || isempty(dataStruct.frameAlignment)
    alignedPoleCoords = [];
    return
end
frameAlignment = dataStruct.frameAlignment;
% if any([getPoleTracks1 getPoleTracks2])
%     if ~isfield(dataStruct,'frameAlignment') || isempty(dataStruct.frameAlignment)
%         alignedPoleCoords = [];
%         return
%     end
%     frameAlignment = dataStruct.frameAlignment;
% end

% get the coords
coordTemp = NaN(dataStruct.dataProperties.movieSize(4),6);
c = 0;
for iPole1 = pole1Idx
    if iPole1 < 0
        c = c + 1;
        feats = pole1Feats(c,:);
        pole1Coords(c).coords = coordTemp;
        frame = 0;
        for iFeat = feats
            frame = frame + 1;
            if iFeat ~= 0
                pole1Coords(c).coords(frame,:) = frameAlignment(frame).alignedCoord(iFeat,:);
            end
        end
    else
        c = c + 1;
        feats = pole1Feats(c,:);
        pole1Coords(c).coords = coordTemp;
        frame = 0;
        for iFeat = feats
            frame = frame + 1;
            if iFeat ~= 0
                pole1Coords(c).coords(frame,:) = frameAlignment(frame).alignedCoord(iFeat,:);
            end
        end
    end
end
c = 0;
for iPole2 = pole2Idx
    if iPole2 < 0
        c = c + 1;
        feats = pole2Feats(c,:);
        pole2Coords(c).coords = coordTemp;
        frame = 0;
        for iFeat = feats
            frame = frame + 1;
            if iFeat ~= 0
                pole2Coords(c).coords(frame,:) = frameAlignment(frame).alignedCoord(iFeat,:);
            end
        end
    else
        c = c + 1;
        feats = pole2Feats(c,:);
        pole2Coords(c).coords = coordTemp;
        frame = 0;
        for iFeat = feats
            frame = frame + 1;
            if iFeat ~= 0
                pole2Coords(c).coords(frame,:) = frameAlignment(frame).alignedCoord(iFeat,:);
            end
        end
    end
end

% save to struct
alignedPoleCoords.pole1Coords = pole1Coords;
alignedPoleCoords.pole2Coords = pole2Coords;


end

