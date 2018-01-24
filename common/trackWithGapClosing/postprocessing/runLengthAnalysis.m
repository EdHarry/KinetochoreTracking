function [runLengthPos,runLengthNeg,dispPos,dispNeg,runLengthPerp,...
    dispPerp,positionProj1,positionProj2] = runLengthAnalysis(tracks,...
    centerCoord,minTrackLen,diffAnalysisRes)
%RUNLENGTHANALYSIS analyzes the speed and time distribution of motion along lineas tracks, relative to a center point
%
%SYNOPSIS [runLengthPos,runLengthNeg,dispPos,dispNeg,runLengthPerp,...
%    dispPerp,positionProj1,positionProj2] = runLengthAnalysis(tracks,...
%    centerCoord,minTrackLen,diffAnalysisRes)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       centerCoord: Coordinates of convergence point of tracks ("center").
%       minTrackLen: Minimum length of a track to be used in analysis.
%                    Optional. Default: 5.
%       diffAnalysisRes: Diffusion analysis results (output of
%                    trackDiffusionAnalysis1). Optional. If not input, will
%                    be calculated here.
%
%OUTPUT runLengthPos: Distribution of run lengths (i.e. number of  
%                     consecutive steps in the same direction) when moving
%                     away from the center. 
%       runLengthNeg: Distribution of run lengths when moving toward the
%                     center.
%       dispPos     : Distribution of displacement magntidues when moving
%                     away from the center.
%       dispNeg     : Distribution of displacement magnitudes when moving
%                     toward the center.
%       runLengthPerp:Distribution of run lengths perpendicular to
%                     direction of motion.
%       dispPerp    : Distribution of displacement magnitudes
%                     perpendicular to direction of motion.
%                     All the above are structure arrays with two entries
%                     (.values), the first for normal diffusion and the
%                     second for super diffusion.
%       positionProj1: Structure array of position projection onto 
%                      preferred direction of motion for linear tracks with
%                      normal diffusion along direction of motion.
%                      For input into armaxFitKalman (can be converted into
%                      input for trajectoryAnalysis using convertTrajectoryData).
%       positionProj2: Structure array of position projection onto 
%                      preferred direction of motion for linear tracks with
%                      super diffusion along direction of motion.
%                      For input into armaxFitKalman (can be converted into
%                      input for trajectoryAnalysis using convertTrajectoryData).
%
%Khuloud Jaqaman, February 2008

%% output
runLengthPos = [];
runLengthNeg = [];
dispPos = [];
dispNeg = [];
runLengthPerp = [];
dispPerp = [];
positionProj1 = [];
positionProj2 = [];

%% input

if nargin < 2 || isempty(tracks) || isempty (centerCoord)
    disp('runLengthAnalysis: Missing input arguments!');
    return
end

if nargin < 3 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 4 || isempty(diffAnalysisRes)
    diffAnalysisRes = trackDiffusionAnalysis1(tracks,1,2,1,[0.05 0.2]);
end

%% preamble

%ignore merges and splits and divide compound tracks back into the
%individual track segments
inputStruct = tracks;
clear tracksFinal
tracks = convStruct2MatIgnoreMS(inputStruct);

%extract track classification from diffAnalysisRes
trackType = vertcat(diffAnalysisRes.classification);

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
tracks = tracks(indx,:);
trackType = trackType(indx,:);

%find linear tracks with normal diffusive behavior along preferred
%direction of motion
tracks1 = tracks(trackType(:,1)==1 & trackType(:,3)==2,:);

%find linear tracks with super-diffusive behavior along preferred direction
%of motion
tracks2 = tracks(trackType(:,1)==1 & trackType(:,3)==3,:);

%find non-linear tracks with normal diffusive behavior
tracks3 = tracks(trackType(:,1)==0 & trackType(:,2)==2,:);

%find non-linear tracks with sub-diffusive behavior
tracks4 = tracks(trackType(:,1)==0 & trackType(:,2)==1,:);

clear criteria tracks indx

%% displacement statistics

%go over both types of linear tracks
for iType = 1 : 4
    
    %get information on current track type and reserve memory
    eval(['tracks = tracks' num2str(iType) ';']);
    numTracks = size(tracks,1);
    positionProjIn = repmat(struct('observations',[],'time',[]),numTracks,1);
    runLengthPosIn = [];
    runLengthNegIn = [];
    dispPosIn = [];
    dispNegIn = [];
    runLengthPerpIn = [];
    dispPerpIn = [];

    %go over all tracks in this category
    for iTrack = 1 : numTracks

        %get the positions in this track and their standard deviations
        %keep NaNs to mark gaps
        trackCoordX = tracks(iTrack,1:8:end)';
        deltaCoordX = tracks(iTrack,5:8:end)';
        trackCoordY = tracks(iTrack,2:8:end)';
        deltaCoordY = tracks(iTrack,6:8:end)';
        trackCoord = [trackCoordX trackCoordY];
        deltaCoord = [deltaCoordX deltaCoordY];

        %project positions onto track's direction of motion
        [posAlongDir,deltaPosAlongDir,velDir] = projectCoordOntoDir(...
            trackCoord,deltaCoord,centerCoord,[]);
        
        %store track projection information in dispProj
        positionProjIn(iTrack).observations = [posAlongDir deltaPosAlongDir];
        positionProjIn(iTrack).time = [(1:length(posAlongDir))' zeros(length(posAlongDir),1)];

        %calculate vector of displacements
        %keep NaNs to mark gaps
        dispVec = trackCoord(2:end,:) - trackCoord(1:end-1,:);

        %calculate the dot product of displacements with the direction
        %vector and its normal
        dispAlongDir = dispVec * velDir;
        dispPerpDir = dispVec * [velDir(2) -velDir(1)]';

        %separate positive displacements along preferred direction from
        %negative displacements along preferred direction
        %NaNs (gaps) have no effect here
        dispPosT = dispAlongDir(dispAlongDir>0);
        dispNegT = -dispAlongDir(dispAlongDir<0);
        
        %collect all displacement magnitudes perpendicular to direction of
        %motion
        dispPerpT = abs(dispPerpDir(~isnan(dispPerpDir)));

        %get the sign of displacements (regardless of value)
        %gaps will constribute NaNs here
        dispSign = sign(dispAlongDir); %along direction
        dispSignPerp = sign(dispPerpDir); %along normal

        %take the difference to find transition points
        dispSignDiff = diff(dispSign); %along direction
        transValue = [0; dispSignDiff(dispSignDiff~=0)];
        transPoints = [0; find(dispSignDiff)];
        dispSignDiffPerp = diff(dispSignPerp); %along normal
        transValuePerp = [0; dispSignDiffPerp(dispSignDiffPerp~=0)];
        transPointsPerp = [0; find(dispSignDiffPerp)];

        %determine run length, i.e. consecutive steps in each direction before
        %switching
        runLength = diff(transPoints); %along preferred direction
        runLengthPerp2 = diff(transPointsPerp); %along normal

        %store run lengths away from nucleus (pos) separate from run lengths
        %toward nucleus (neg)
        %remove run lengths involving gaps (NaNs)
        runLengthPosT = []; %along preferred direction
        runLengthNegT = [];
        for iTrans = 1 : length(runLength)
            if ~isnan(transValue(iTrans+1)) && ~isnan(transValue(iTrans))
                if transValue(iTrans+1) > 0
                    runLengthNegT = [runLengthNegT; runLength(iTrans)];
                else
                    runLengthPosT = [runLengthPosT; runLength(iTrans)];
                end
            end
        end
        runLengthPerpT = []; %along normal
        for iTrans = 1 : length(runLengthPerp2)
            if ~isnan(transValuePerp(iTrans+1)) && ~isnan(transValuePerp(iTrans))
                runLengthPerpT = [runLengthPerpT; runLengthPerp2(iTrans)];
            end
        end

        %separate positive steps from negative steps
        %add this track's results to the rest
        runLengthPosIn = [runLengthPosIn; runLengthPosT];
        runLengthNegIn = [runLengthNegIn; runLengthNegT];
        dispPosIn = [dispPosIn; dispPosT];
        dispNegIn = [dispNegIn; dispNegT];
        runLengthPerpIn = [runLengthPerpIn; runLengthPerpT];
        dispPerpIn = [dispPerpIn; dispPerpT];

    end %(for iTrack = 1 : numTracks)
    
    %store results for output
    runLengthPos(iType).values = runLengthPosIn;
    runLengthNeg(iType).values = runLengthNegIn;
    dispPos(iType).values = dispPosIn;
    dispNeg(iType).values = dispNegIn;
    runLengthPerp(iType).values = runLengthPerpIn;
    dispPerp(iType).values = dispPerpIn;
    eval(['positionProj' num2str(iType) ' = positionProjIn;'])
    
end %(for iType = 1 : 2)

%% ~~~ the end ~~~



