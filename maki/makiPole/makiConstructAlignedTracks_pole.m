function alignedTracks = makiConstructAlignedTracks_pole( dataStruct )
%MAKICONSTRUCTALIGNEDTRACKS Gives Aligned Coordinates of all Tracks
%alignedTracks - track structure with the x,y,z data taken from the
%alignedCoordinates
%
%atLeastOnePlaneFit - flag for whether or not there was at least one
%sucessfull plane fit in the movie
% EHarry Jan 2012

%planeFit = dataStruct.planeFit; % copy out the nessesary fields
frameAlignment = dataStruct.frameAlignment_pole;
tracks = dataStruct.tracks_pole;

if isempty(frameAlignment) || isempty(tracks) % end if data is missing
    alignedTracks = [];
    return
end

alignedTracks = tracks; % copy the tracks structure

% if isempty(cat(1,planeFit.planeVectors)) % look for at least one planeFit by looking for at least one non-empty list of planeVectors
%     atLeastOnePlaneFit = 0;
% else
%     atLeastOnePlaneFit = 1;
% end

for iTrack = 1:length(tracks)
    
    track = tracks(iTrack);
    
    featIdx = track.tracksFeatIndxCG; % featureIdx of the track
    
    trackStart = track.seqOfEvents(1,1); % start frame of the track
    trackEnd = track.seqOfEvents(2,1); % end frame of the track
    
    featCount = 0; % feature idx counter
    
    for iFrame = trackStart:trackEnd % loop over frames where the track exists
        
        featCount = featCount + 1; % update the feat counter
        
        feat = featIdx(featCount); % get the feature idx
        
        if feat ~= 0 % if the feature exists
            
            xPos = (featCount - 1).*8 + 1;
            yPos = (featCount - 1).*8 + 2;
            zPos = (featCount - 1).*8 + 3;
            
            dxPos = (featCount - 1).*8 + 5;
            dyPos = (featCount - 1).*8 + 6;
            dzPos = (featCount - 1).*8 + 7;
            
            alignedTracks(iTrack).tracksCoordAmpCG(xPos) = frameAlignment(iFrame).alignedCoord(feat,1); % get the data from the alignedCoords
            alignedTracks(iTrack).tracksCoordAmpCG(yPos) = frameAlignment(iFrame).alignedCoord(feat,2);
            alignedTracks(iTrack).tracksCoordAmpCG(zPos) = frameAlignment(iFrame).alignedCoord(feat,3);
            
            alignedTracks(iTrack).tracksCoordAmpCG(dxPos) = frameAlignment(iFrame).alignedCoord(feat,4);
            alignedTracks(iTrack).tracksCoordAmpCG(dyPos) = frameAlignment(iFrame).alignedCoord(feat,5);
            alignedTracks(iTrack).tracksCoordAmpCG(dzPos) = frameAlignment(iFrame).alignedCoord(feat,6);
            
        end
        
    end
    
end

end

