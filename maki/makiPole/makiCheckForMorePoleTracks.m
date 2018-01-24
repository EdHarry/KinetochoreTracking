function [newPoleTracks,found] = makiCheckForMorePoleTracks( dataStruct )
%MAKICHECKFORMOREPOLETRACKS checks pole_tracks for the possibility of finding more spindle pole tracks
% EHarry Jan 2012

found = 0; % indicator of whether or not extra poles have been found

minTrackLengthPercent = 50; % minimum percentage length of a track w.r.t the movie length to be considred for a pole

noTimePoints = dataStruct.dataProperties.movieSize(4); % number of time points


if isfield(dataStruct,'poles')
    poles = dataStruct.poles; % copy out the tracks
else
    newPoleTracks=[];
    return
end


if isempty(poles)
    newPoleTracks=[];
    return
end



if isfield(dataStruct,'tracks_pole')
    tracks_pole = dataStruct.tracks_pole; % copy out the tracks
else
    newPoleTracks=[];
    return
end


if isempty(tracks_pole)
    newPoleTracks=[];
    return
end


trackStats = catStruct(3,'tracks_pole.seqOfEvents');
trackLength = squeeze(trackStats(2,1,:)-trackStats(1,1,:)+1);

goodTracks = find(trackLength >= (minTrackLengthPercent./100).*noTimePoints);
nGoodTracks = length(goodTracks);



maxCentriolAvDis = 1; % max av. distance a pole's other centriol can be away from itself
minCentriolCrossCorr = 0.9; % minimum crossCorr of two centriols at lag 0

costCentriol = -ones(2,nGoodTracks); % cost matix for centriol id

for i = 1:length(poles.pole1Track)
    pole1Complete_temp(:,:,i) = getTrackCoords(poles.pole1Track(i) , noTimePoints);
end

for i = 1:length(poles.pole2Track)
    pole2Complete_temp(:,:,i) = getTrackCoords(poles.pole2Track(i) , noTimePoints);
end

trackCoords_1(:,1:3) = nanmean(pole1Complete_temp(:,1:3,:),3);
trackCoords_2(:,1:3) = nanmean(pole2Complete_temp(:,1:3,:),3);


%trackCoords_1 = getTrackCoords(poles.pole1Track , noTimePoints); % pole properties
track1Norm = normList(trackCoords_1);
%trackCoords_2 = getTrackCoords(poles.pole2Track , noTimePoints);
track2Norm = normList(trackCoords_2);

if length(poles.pole1TrackIdx) == 1 || length(poles.pole2TrackIdx) == 1 % only if some centrols haven't been found
    for i = 1:nGoodTracks
        
        
        trackCoords_i = getTrackCoords(tracks_pole(goodTracks(i)) , noTimePoints);
        trackNorm_i = normList(trackCoords_i);
        
        distances_1 = getDistances(trackCoords_1,trackCoords_i);% get the distances between the tracks
        distances_2 = getDistances(trackCoords_2,trackCoords_i);% get the distances between the tracks
        
        avDis_1 = nanmean(distances_1); % get the averages, ignoring NaNs
        avDis_2 = nanmean(distances_2);
        
        gamma_1 = crossCorr(track1Norm,trackNorm_i,0); % crossCorr between the tracks at lag 0
        gamma_1 = gamma_1(1);
        
        gamma_2 = crossCorr(track2Norm,trackNorm_i,0); % crossCorr between the tracks at lag 0
        gamma_2 = gamma_2(1);
        
        % for pole1
        if avDis_1 < maxCentriolAvDis && gamma_1 > minCentriolCrossCorr && length(poles.pole1TrackIdx) == 1
            costCentriol(1,i) = 1./gamma_1; % add to cost matrix, cost is the inverse of the cross corr
        end
        
        % for pole2
        if avDis_2 < maxCentriolAvDis && gamma_2 > minCentriolCrossCorr && length(poles.pole2TrackIdx) == 1
            costCentriol(2,i) = 1./gamma_2; % add to cost matrix, cost is the inverse of the cross corr
        end
        
    end
end

if ~all(costCentriol(:)==-1) % if links are possible
    %     dataStruct.poles = poles;
    %     return
    %end
    
    centriolLinks = lap(costCentriol,[],[],1);
    centriolLinks = centriolLinks(1:2);
    
    % Get Centriols if Links are Good
    if centriolLinks(1) <= nGoodTracks
        % remove dependencies structure if it's part of the tracks
        % structure and not part of poles structure
        pole1Dep = isfield(poles.pole1Track, 'dependencies');
        tracks_poleDep = isfield(tracks_pole, 'dependencies');
        
        if pole1Dep ~= tracks_poleDep
            if pole1Dep
                poles.pole1Track = rmfield(poles.pole1Track, 'dependencies');
            else
                tracks_pole = rmfield(tracks_pole, 'dependencies');
            end
        end
        
        poles.pole1Track(2) = tracks_pole(goodTracks(centriolLinks(1)));
        poles.pole1TrackIdx(2) = -goodTracks(centriolLinks(1)); % minus indicates that the track is from tracks_pole
        found = found + 1;
    end
    
    if centriolLinks(2) <= nGoodTracks
        pole2Dep = isfield(poles.pole2Track, 'dependencies');
        tracks_poleDep = isfield(tracks_pole, 'dependencies');
        
        if pole2Dep ~= tracks_poleDep
            if pole2Dep
                poles.pole2Track = rmfield(poles.pole2Track, 'dependencies');
            else
                tracks_pole = rmfield(tracks_pole, 'dependencies');
            end
        end
        
        
        poles.pole2Track(2) = tracks_pole(goodTracks(centriolLinks(2)));
        poles.pole2TrackIdx(2) = -goodTracks(centriolLinks(2));
        found = found + 1;
    end
end

if found > 0
    found = 1;
end

newPoleTracks = poles;








%%  SUBFUNTIONS

    function trackCoords = getTrackCoords(track , noTimePoints)
        % takes a track struct from maki and returns an nx3 matirx of the
        % track coords (absolute coords)
        
        trackCoords = NaN(noTimePoints,3); % predefine the coords, NaNs will be left for missing time points
        
        startTime = track.seqOfEvents(1,1); % get track start and end times
        endTime = track.seqOfEvents(2,1);
        
        x = track.tracksCoordAmpCG(1:8:end); % get the coords
        y = track.tracksCoordAmpCG(2:8:end);
        z = track.tracksCoordAmpCG(3:8:end);
        
        trackCoords(startTime:endTime,:) = [x;y;z]';
    end

    function distances = getDistances(trackCoords1,trackCoords2)
        % gets the distances between two tracks
        
        distances = normList(trackCoords1 - trackCoords2); % ditances are just the normals of the difference between the tracks
    end

end

