function dataStruct = makiGenerateTracks_pole(dataStruct)
% EHarry Jan 2012

%% ORIGINAL HEADER
% % %MAKIGENERATETRACKS tracks kinetochores throughout a movie
% % %
% % %SYNOPSIS dataStruct = makiGenerateTracks(dataStruct)
% % %
% % %INPUT  dataStruct: dataStruct as in makiMakeDataStruct with the
% % %                   fields "dataProperties", "initCoord" & "planeFit_pole".
% % %                   Field "planeFit_pole" can be empty.
% % %                   Optional. Loaded interactively if not input.
% % %
% % %OUTPUT dataStruct: Same as input, with added field "tracks"
% % %
% % %Khuloud Jaqaman, July 2007
% % %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load dataStruct if not input
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end

if ~isfield(dataStruct,'poleReferenceFrame') || isempty(dataStruct.poleReferenceFrame);
    dataStruct.tracks_pole = [];
    return
end


new2OriginalIdxMap = []; % this will serve as a flag for using polar coords or not

%get number of time points in movie
nTimepoints = dataStruct.dataProperties.movieSize(4);

%get kinetochore coordinates and amplitude
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[]),...
    nTimepoints,1);

% if polar coordinates are to be used ...
if dataStruct.dataProperties.tracksParam.pole ==1 && ~isempty(dataStruct.poleReferenceFrame)
    
    clear new2OriginalIdxMap % clear the empty vector
    new2OriginalIdxMap(1:nTimepoints) = struct('map',[]);% index map back to original as is in initCoord
    
    %get the polar coordinates
    for iTime = 1 : nTimepoints
        allCoord = dataStruct.poleReferenceFrame(iTime).poleCoords_cartisian;
        goodCoordIdx = ~isnan(allCoord(:,1)); % index to goodCoords (not poles or frames with no poles)
        new2OriginalIdxMap(iTime).map = find(goodCoordIdx); 
        allCoord = allCoord(goodCoordIdx,:); % new coords with no NaNs
        movieInfo(iTime).xCoord = [allCoord(:,1) allCoord(:,4)];
        movieInfo(iTime).yCoord = [allCoord(:,2) allCoord(:,5)];
        movieInfo(iTime).zCoord = [allCoord(:,3) allCoord(:,6)];
        amp = dataStruct.initCoord(iTime).amp;
        movieInfo(iTime).amp = amp(goodCoordIdx,:);
    end
    
    
    %if rotated coordinates are to be used ...    
elseif dataStruct.dataProperties.tracksParam.rotate == 1 && ~isempty(dataStruct.planeFit_pole)
    
    %get the rotated coordinates
    for iTime = 1 : nTimepoints
        allCoord = dataStruct.planeFit_pole(iTime).rotatedCoord;
        movieInfo(iTime).xCoord = [allCoord(:,1) allCoord(:,4)];
        movieInfo(iTime).yCoord = [allCoord(:,2) allCoord(:,5)];
        movieInfo(iTime).zCoord = [allCoord(:,3) allCoord(:,6)];
        movieInfo(iTime).amp = dataStruct.initCoord(iTime).amp;
    end
    
else %if the original coordinates are to be used
    
    %get the original coordinates
    for iTime = 1 : nTimepoints
        allCoord = dataStruct.initCoord(iTime).allCoord;
        movieInfo(iTime).xCoord = [allCoord(:,1) allCoord(:,4)];
        movieInfo(iTime).yCoord = [allCoord(:,2) allCoord(:,5)];
        movieInfo(iTime).zCoord = [allCoord(:,3) allCoord(:,6)];
        movieInfo(iTime).amp = dataStruct.initCoord(iTime).amp;
    end
    
    %calculate the center of mass in each frame
    centerOfMass = zeros(nTimepoints,3);
    for iTime = 1 : nTimepoints
       centerOfMass(iTime,:) = [mean(movieInfo(iTime).xCoord(:,1)) ...
           mean(movieInfo(iTime).yCoord(:,1)) mean(movieInfo(iTime).zCoord(:,1))];
    end
    
    %shift coordinates by center of mass to make the origin in each frame 
    %at its center of mass
    for iTime = 1 : nTimepoints
        movieInfo(iTime).xCoord(:,1) = movieInfo(iTime).xCoord(:,1) - centerOfMass(iTime,1);
        movieInfo(iTime).yCoord(:,1) = movieInfo(iTime).yCoord(:,1) - centerOfMass(iTime,2);
        movieInfo(iTime).zCoord(:,1) = movieInfo(iTime).zCoord(:,1) - centerOfMass(iTime,3);
    end
    
end
    
%get number of features in each frame
if ~isfield(movieInfo,'num')
    for iTime = 1 : nTimepoints
        movieInfo(iTime).num = size(movieInfo(iTime).xCoord,1);
    end
end

%collect coordinates and their std in one matrix in each frame
if ~isfield(movieInfo,'allCoord')
    for iTime = 1 : nTimepoints
        movieInfo(iTime).allCoord = [movieInfo(iTime).xCoord ...
            movieInfo(iTime).yCoord movieInfo(iTime).zCoord];
    end
end

%calculate nearest neighbor distance for each feature in each frame
if ~isfield(movieInfo,'nnDist')

    for iTime = 1 : nTimepoints
        
        switch movieInfo(iTime).num

            case 0 %if there are no features

                %there are no nearest neighbor distances
                nnDist = zeros(0,1);

            case 1 %if there is only 1 feature

                %assign nearest neighbor distance as 1000 pixels (a very big
                %number)
                nnDist = 1000;

            otherwise %if there is more than 1 feature

                %compute distance matrix
                nnDist = createDistanceMatrix(movieInfo(iTime).allCoord(:,1:2:end),...
                    movieInfo(iTime).allCoord(:,1:2:end));

                %sort distance matrix and find nearest neighbor distance
                nnDist = sort(nnDist,2);
                nnDist = nnDist(:,2);

        end

        %store nearest neighbor distance
        movieInfo(iTime).nnDist = nnDist;

    end
    
end

%get kinetochore classification in each frame if available
if ~isempty(dataStruct.planeFit_pole) %if the plane fit has been done
    
    %extract planeFit_pole field from structure
    planeFit_pole = dataStruct.planeFit_pole;
    
    %assign kinetochore types per frame
    for iTime = 1 : nTimepoints
        
        un = planeFit_pole(iTime).unalignedIdx;
        lag = planeFit_pole(iTime).laggingIdx;
        
        if ~isempty(new2OriginalIdxMap) % if polar coord have ben used have to convert idexes to new system
            un = un(ismember(un,new2OriginalIdxMap(iTime).map));
            lag = lag(ismember(lag,new2OriginalIdxMap(iTime).map));
            un = find(ismember(new2OriginalIdxMap(iTime).map,un));
            lag = find(ismember(new2OriginalIdxMap(iTime).map,lag));
        end
        
        
        kinType = zeros(movieInfo(iTime).num,1); %inlier
        kinType(un) = 1; %unaligned
        kinType(lag) = 2; %lagging
        movieInfo(iTime).kinType = kinType;
    end
    
else %if not
    
    %treat all kinetochores as inliers
    for iTime = 1 : nTimepoints
        movieInfo(iTime).kinType = zeros(movieInfo(iTime).num,1);
    end

end

%get tracking parameters
gapCloseParam = dataStruct.dataProperties.tracksParam.gapCloseParam;
costMatrices = dataStruct.dataProperties.tracksParam.costMatrices;
kalmanFunctions = dataStruct.dataProperties.tracksParam.kalmanFunctions;

%call tracker
try
    
    %track the kinetochores
    tracks = trackCloseGapsKalman(movieInfo,...
        costMatrices,gapCloseParam,kalmanFunctions,3,0,0);

    %replace the coordinate used for tracking (whether the rotated
    %coordinates or the original coordinates shifted by center of mass) by
    %the original coordinates
    for iTrack = 1 : length(tracks)

        %store coordinates used for tracking in another field
        tracks(iTrack).coordAmp4Tracking = tracks(iTrack).tracksCoordAmpCG;

        %fetch the start and end time of this track
        startTime = tracks(iTrack).seqOfEvents(1,1);
        endTime = tracks(iTrack).seqOfEvents(2,1);
        
        %go over all frames where this track exists
        for iFrame =  startTime : endTime
            
            
            %get the feature making up this track in this frame
            iFeature = tracks(iTrack).tracksFeatIndxCG(iFrame-startTime+1);
            
            
            %if there is a feature (not a gap)
            if iFeature ~= 0
                
                if ~isempty(new2OriginalIdxMap) % if polar coord have been used then iFeature needs to be mapped back to the originals
                    iFeature = new2OriginalIdxMap(iFrame).map(iFeature);
                    tracks(iTrack).tracksFeatIndxCG(iFrame-startTime+1) = iFeature; % write the original feature index to the track
                end
                
                
                %replace coordiantes and their stds
                tracks(iTrack).tracksCoordAmpCG(1,(iFrame-startTime)*8+1:...
                    (iFrame-startTime)*8+3) = dataStruct.initCoord(iFrame).allCoord(iFeature,1:3);
                tracks(iTrack).tracksCoordAmpCG(1,(iFrame-startTime)*8+5:...
                    (iFrame-startTime)*8+7) = dataStruct.initCoord(iFrame).allCoord(iFeature,4:6);

            end

        end

    end %(for iTrack = 1 : numTracks)

catch
    disp('error in tracking')
end



%% dependencies
% tracks_pole depends on dataProperites, initCoords, planeFit, tracks,
% sisterList, updatedClass, frameAlignment, poles, poleReferenceFrame and planeFit_pole
% dependencies = struct('dataProperties',[],'initCoord',[],'planeFit',[],'tracks',[],'sisterList',[],'updatedClass',[],'frameAlignment',[],'poles',[],'poleReferenceFrame',[],'planeFit_pole',[]);

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
polesName = dataStruct.polesName;
polesV = getVersion(polesName);
poleReferenceFrameName = dataStruct.poleReferenceFrameName;
poleReferenceFrameV = getVersion(poleReferenceFrameName);
planeFit_poleName = dataStruct.planeFit_poleName;
planeFit_poleV = getVersion(planeFit_poleName);


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
dependencies.poles = polesV;
dependencies.poleReferenceFrame = poleReferenceFrameV;
dependencies.planeFit_pole = planeFit_poleV;

% save into the first tracks_pole
tracks(1).dependencies = dependencies;




%store tracks in dataStruct
dataStruct.tracks_pole = tracks;
