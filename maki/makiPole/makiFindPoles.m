function dataStruct = makiFindPoles( dataStruct )
%MAKIFINDPOLES Findes pole tracks in a maki movie
%   Input:  dataStruct:     maki data-structure with fields: tracks and
%                           possibly planeFit
%
%   Output: dataStruct:     maki data-structure with additional field:
%   poles
%
% eharry, 02/12/11


%%%%  MAIN

%% Initial Stuff

minTrackLengthPercent = 50; % minimum percentage length of a track w.r.t the movie length to be considred for a pole
avarageLengthOfSpindle = 11; % value for the average length of the spindle
minAvLength = 8; % minimum av distance between pole tracks
maxAvAngle = 0.698131700797732;  % maximum av angle that the pole track axis makes with the planeFit, if it exists (40 deg in rad)

noTimePoints = dataStruct.dataProperties.movieSize(4); % number of time points


if isfield(dataStruct,'tracks')
    tracks = dataStruct.tracks; % copy out the tracks
else
    dataStruct.poles=[];
    return
end

if isfield(dataStruct,'frameAlignment')
    frameAlignment = dataStruct.frameAlignment; % copy out frameAlignment if it exists
else
    dataStruct.poles=[];
    return
end

if isfield(dataStruct,'updatedClass')
    updatedClass = dataStruct.updatedClass; % copy out frameAlignment if it exists
else
    dataStruct.poles=[];
    return
end

if isempty(tracks)
    dataStruct.poles=[];
    return
end

if isempty(frameAlignment)
    dataStruct.poles=[];
    return
end

if isempty(updatedClass)
    dataStruct.poles=[];
    return
end

trackStats = catStruct(3,'tracks.seqOfEvents');
trackLength = squeeze(trackStats(2,1,:)-trackStats(1,1,:)+1);

goodTracks = find(trackLength >= (minTrackLengthPercent./100).*noTimePoints);
nGoodTracks = length(goodTracks);

tracks_original = tracks; % save the tracks as an original

[ tracks , atLeastOnePlaneFit ] = makiConstructAlignedTracks( dataStruct ); % get the aligned tracks

phase = catStruct(1,'updatedClass.phase',[],1); % phases of the cell cycle

firstFrameAPhase = find(phase=='a',1); % find the first frame theat the cell goes into aphase, if it does at all

if isempty(firstFrameAPhase) % if the cell did not go into aphase
    endTimeToTest = noTimePoints; % test tracks until the end of the movie
else
    endTimeToTest = firstFrameAPhase - 1; % else test upto the time just before aphase
end

if atLeastOnePlaneFit % if there was at least one planeFit then the normal angles can all be considerd to be [1 0 0]
    normals = repmat([1 0 0],noTimePoints,1);
else % otherwise set the normals to NaN (do not test the angles)
    normals = NaN(noTimePoints,3);
end

% if ~isempty(frameAlignment)
%     normals = getNormals(frameAlignment); % get the normls of the metaphse plate, if it exists
% else
%     normals = NaN(length(frameAlignment),3);
% end

%% Main Cost Matrix Loop

cost = NaN(nGoodTracks); % initilise cost matrix for pairing

for i = 1:nGoodTracks % loop over all pairs of tracks
    for j = i+1:nGoodTracks
        
        %display([int2str(goodTracks(i)) ' -> ' int2str(goodTracks(j))]);
        
        trackCoords_i = getTrackCoords(tracks(goodTracks(i)) , noTimePoints); % get the trackCoords
        trackCoords_j = getTrackCoords(tracks(goodTracks(j)) , noTimePoints);
        
        featIdx_i = getFeatIdx(tracks(goodTracks(i)),noTimePoints);
        featIdx_j = getFeatIdx(tracks(goodTracks(j)),noTimePoints);
        
        distances = getDistances(trackCoords_i,trackCoords_j);% get the distances between the tracks
        angles = getAngleWithNormal(trackCoords_i,trackCoords_j,normals);% get the angles with the normal, if it exists
        
        %         avDis = nanmean(distances); % get the averages, ignoring NaNs
        %         avAngle = nanmean(angles);
        
        % make this robust
        avDis = robustMean(distances); % get the averages, ignoring NaNs
        avAngle = robustMean(angles);
        
        %display(['avDis_' int2str(goodTracks(i)) '->' int2str(goodTracks(j)) ' = ' num2str(avDis)]);
        %display(['avAngle_' int2str(goodTracks(i)) '->' int2str(goodTracks(j)) ' = ' num2str(avAngle)]);
        
        if avDis < minAvLength || avAngle > maxAvAngle
            % do nothing
        else
            % potentially add to cost matrix
            
            toAdd = 1;
            isUnaligned = 1; % flag for both tracks being unaligned at all timepoints with a plate fit
            position_i = NaN; % marker for the side of the metaphse plate (if) it exists for track_i
            position_j = NaN; % marker for the side of the metaphse plate (if) it exists for track_j
            
            for t = 1:endTimeToTest
                
                if ~isnan(normals(t,1)) % if a plate exist at this time point
                    
                    if ~isnan(featIdx_i(t))  % if track_i exists at this point
                        isUnaligned = isUnaligned && ismember(featIdx_i(t),updatedClass(t).unalignedIdx);
                        
                        if isnan(position_i)
                            position_i = sign(frameAlignment(t).alignedCoord(featIdx_i(t),1)); % only record the first available side of the plate for each track, a track could cross the plate, but then it would become aligned and would be pickup up by the other logic check
                        end
                    end
                    
                    if ~isnan(featIdx_j(t))  % if track_j exists at this point
                        isUnaligned = isUnaligned && ismember(featIdx_j(t),updatedClass(t).unalignedIdx);
                        
                        
                        if isnan(position_j)
                            position_j = sign(frameAlignment(t).alignedCoord(featIdx_j(t),1));
                        end
                    end
                    
                end
                
                if ~isUnaligned || position_i == position_j
                    toAdd = NaN; % NaN if one of the tracks is not unaligned or if both tracks are on the same side of the plate, if it exists
                    break
                end
                
            end
            
            if isnan(avAngle)
                avAngle = 1; % set the av angle to 1 if no angles existed (no planeFit)
            end
            
            
            cost(i,j) = abs(avDis - avarageLengthOfSpindle).*avAngle.*toAdd; % cost is the product of av dis minus the av length of the spindle, the avAngle with the normal, and the flag of whether a track was unaligned or not
            cost(j,i) = cost(i,j); % cost is symmetric
            
        end
        
    end
end



cost(isnan(cost)) = -1; % set NaNs i.e. impossible links to -1 for the lap


if all(all(cost==-1)) % if no poles possible
    dataStruct.poles = [];
    return
end

%% LAP

[r2c,c2r] = lap(cost,[],[],1); % try the lap


r2c = double(r2c(1:nGoodTracks));
r2c(r2c>nGoodTracks) = NaN;
c2r = double(c2r(1:nGoodTracks));
c2r(c2r>nGoodTracks) = NaN;

goodPairIdxL = r2c==c2r;



%% RESOLVE POLYGONS

% identify polygons. Polygons have to be closed, thus, it should not matter
% where we start. Also, since the distance between polygons and the rest is
% hopefully fairly large, we don't care about neighborhood.
polygonIdx = find(~goodPairIdxL);

% remove the not-linked tracks
polygonIdx(isnan(r2c(polygonIdx))) = [];
polyList = [];
while ~isempty(polygonIdx)
    polyList(1) = polygonIdx(1);
    polygonIdx(1) = [];
    done = false;
    while ~done
        % look up the row the last corner links to
        nextCorner = r2c(polyList(end));
        % check whether the new corner has already been used
        if any(nextCorner == polyList)
            % if yes, exit. The polygon is complete
            done = true;
        else
            polyList(end+1) = nextCorner; %#ok<AGROW>
            % remove corner from polygonIdx
            polygonIdx(polygonIdx==nextCorner) = [];
        end
    end % identify polygon
    
    % within the polygon: find closest distance to identify first pair.
    % Remove it, and check for more pairs. This will potentially result in
    % more pairs than removal of large distances starting from a tetragon
    done = false;
    while ~done
        % read current cost matrix
        currentCost = cost(polyList,polyList);
        currentCost(currentCost==-1) = inf;
        % find pair with lowest cost
        [~,minIdx] = min(currentCost(:));
        [idx1,idx2] = ind2sub(size(currentCost),minIdx);
        % write pair into r2c
        r2c(polyList(idx1)) = polyList(idx2);
        r2c(polyList(idx2)) = polyList(idx1);
        c2r(polyList(idx1)) = polyList(idx2);
        c2r(polyList(idx2)) = polyList(idx1);
        
        % check whether there are still tracks to link
        polyList([idx1,idx2]) = [];
        
        if length(polyList) > 1 && any(~isinf(currentCost(:)))
            % continue
        else
            % clear polyList, write NaN into r2c, c2r
            r2c(polyList) = NaN;
            c2r(polyList) = NaN;
            polyList = [];
            done = true;
        end
    end % resolve individual polygons
    
end % resolve all polygons



goodPairIdxL = r2c==c2r;
linkedIdx=sub2ind([nGoodTracks nGoodTracks],find(goodPairIdxL),r2c(goodPairIdxL));
% trackPairs = [goodTracks(goodPairIdxL),goodTracks(r2c(goodPairIdxL))];
% trackPairs = sort(trackPairs,2);
% trackPairs = unique(trackPairs,'rows');
trackPairs = [goodTracks(goodPairIdxL),goodTracks(r2c(goodPairIdxL)),cost(linkedIdx)];

% remove redundancy
trackPairs(:,1:2) = sort(trackPairs(:,1:2),2);
trackPairs = unique(trackPairs,'rows');

%% Identify Poles and Centriols if Possible

%noPairs = size(trackPairs,1);

% if noPairs == 1 % only one pair of tracks identified
%
%     poles.otherCentroil = [0 0]; % indicates that neither pole has its other centriol
%
%     poles.pole1Track = tracks(trackPairs(1));
%     poles.pole2Track = tracks(trackPairs(2));
%     poles.pole1TrackIdx = trackPairs(1);
%     poles.pole2TrackIdx = trackPairs(2);
%
%     dataStruct.poles = poles;
%
%     return
%
% else % more than one pair to consider, possiblly have other centriols
%
%
%
% end

[~,polesIdx] = min(trackPairs(:,3)); % take the pair with the lowest cost as the poles base
pole1Idx = trackPairs(polesIdx,1);
pole2Idx = trackPairs(polesIdx,2);

poles.pole1Track = tracks(pole1Idx);  % store initial information
poles.pole2Track = tracks(pole2Idx);
poles.pole1TrackIdx = pole1Idx;
poles.pole2TrackIdx = pole2Idx;

maxCentriolAvDis = 1; % max av. distance a pole's other centriol can be away from itself
minCentriolCrossCorr = 0.9; % minimum crossCorr of two centriols at lag 0

costCentriol = NaN(2,nGoodTracks); % cost matix for centriol id

trackCoords_1 = getTrackCoords(tracks(pole1Idx) , noTimePoints); % pole properties
track1Norm = normList(trackCoords_1);
trackCoords_2 = getTrackCoords(tracks(pole2Idx) , noTimePoints);
track2Norm = normList(trackCoords_2);

for i = 1:nGoodTracks
    if goodTracks(i) ~= pole1Idx && goodTracks(i) ~= pole2Idx
        
        trackCoords_i = getTrackCoords(tracks(goodTracks(i)) , noTimePoints);
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
        if avDis_1 < maxCentriolAvDis && gamma_1 > minCentriolCrossCorr
            costCentriol(1,i) = 1./gamma_1; % add to cost matrix, cost is the inverse of the cross corr
        end
        
        % for pole2
        if avDis_2 < maxCentriolAvDis && gamma_2 > minCentriolCrossCorr
            costCentriol(2,i) = 1./gamma_2; % add to cost matrix, cost is the inverse of the cross corr
        end
    end
end

%% LAP on Centriols

costCentriol(isnan(costCentriol)) = -1;

if ~all(costCentriol(:)==-1) % if links are possible
    %     dataStruct.poles = poles;
    %     return
    %end
    
    centriolLinks = lap(costCentriol,[],[],1);
    centriolLinks = centriolLinks(1:2);
    
    % Get Centriols if Links are Good
    if centriolLinks(1) <= nGoodTracks
        poles.pole1Track(2) = tracks(goodTracks(centriolLinks(1)));
        poles.pole1TrackIdx(2) = goodTracks(centriolLinks(1));
    end
    
    if centriolLinks(2) <= nGoodTracks
        poles.pole2Track(2) = tracks(goodTracks(centriolLinks(2)));
        poles.pole2TrackIdx(2) = goodTracks(centriolLinks(2));
    end
    
end

% put pole1 on the +ve side of the plate

if poles.pole1Track(1).tracksCoordAmpCG(1) < 0
    temp.track = poles.pole1Track;
    temp.idx = poles.pole1TrackIdx;
    poles.pole1Track = poles.pole2Track;
    poles.pole1TrackIdx = poles.pole2TrackIdx;
    poles.pole2Track = temp.track;
    poles.pole2TrackIdx = temp.idx;
end


% put back the original track info

for i = 1:length(poles.pole1TrackIdx)
    poles.pole1Track(i) = tracks_original(poles.pole1TrackIdx(i));
end

for i = 1:length(poles.pole2TrackIdx)
    poles.pole2Track(i) = tracks_original(poles.pole2TrackIdx(i));
end


%% dependencies
% poles depends on dataProperites, initCoords, planeFit, tracks,
% sisterList, updatedClass and frameAlignment
% dependencies = struct('dataProperties',[],'initCoord',[],'planeFit',[],'tracks',[],'sisterList',[],'updatedClass',[],'frameAlignment',[]);

dataPropName = dataStruct.dataPropertiesName;
dataPropV = getVersion(dataPropName);
initCoordName = dataStruct.initCoordName;
initCoordV = getVersion(initCoordName);
planeFitName = dataStruct.planeFitName;
planeFitV = getVersion(planeFitName);
tracksName = dataStruct.tracksName;
tracksV = getVersion(tracksName);
sisterListName = dataStruct.sisterListName;
sisterListV = getVersion(sisterListName);
updatedClassName = dataStruct.updatedClassName;
updatedClassV = getVersion(updatedClassName);
frameAlignmentName = dataStruct.frameAlignmentName;
frameAlignmentV = getVersion(frameAlignmentName);

dependencies.dataProperties = dataPropV;
dependencies.initCoord = initCoordV;
if ~isempty(dataStruct.planeFit)
    dependencies.planeFit = planeFitV;
end
if ~isempty(dataStruct.tracks)
    dependencies.tracks = tracksV;
end
if ~isempty(dataStruct.sisterList)
    dependencies.sisterList = sisterListV;
end
if ~isempty(dataStruct.updatedClass)
    dependencies.updatedClass = updatedClassV;
end
if ~isempty(dataStruct.frameAlignment)
    dependencies.frameAlignment = frameAlignmentV;
end

% save into the first poles
poles(1).dependencies = dependencies;


dataStruct.poles = poles;

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


    function angles = getAngleWithNormal(trackCoords1,trackCoords2,normals)
        % takes the coords of two tracks and gives the angles (radians) of the vector
        % between then and the normal vectors
        
        trackVectors = trackCoords1 - trackCoords2; % take the vecotrs between the tracks
        
        [~,normals] = normList(normals); % normalise the normals
        [~,normedVectors]=normList(trackVectors); % normalise the track vecotrs
        
        angles = acos(dot(normedVectors,normals,2)); % calcule the angles
        
        angles(angles>pi/2) = pi - angles(angles>pi/2); % correct for incorrect direction of vectors
    end

%     function normals = getNormals(planeFit)
%         % gets the normals to the planeFit as a nx3 matrix, returns NaNs
%         % for times with no normals
%
%         normals = NaN(length(planeFit),3);
%
%         for ii = 1:length(planeFit)
%             vectors = planeFit(ii).planeVectors;
%             if ~isempty(vectors)
%                 normals(ii,:) = vectors(:,1)';
%             end
%         end
%     end

    function distances = getDistances(trackCoords1,trackCoords2)
        % gets the distances between two tracks
        
        distances = normList(trackCoords1 - trackCoords2); % ditances are just the normals of the difference between the tracks
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

