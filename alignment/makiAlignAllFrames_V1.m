function dataStruct = makiAlignAllFrames_V1(dataStruct)
% edit of makiAlignFrames to align all frames to the frame with the "best"
% plane-fit
% EHarry Nov 2012
%% original header
%MAKIROTATEFRAMES aligns frames based on translation and rotation estimates derived from tracks
%
%SYNOPSIS dataStruct = makiAlignFrames(dataStruct)
%
%INPUT  dataStruct : dataStruct as in makiMakeDataStruct with the
%                    fields "dataProperties", "initCoord", "planeFit",
%                    "tracks" & "updatedClass". Field "planeFit" can be
%                    empty.
%                    Optional. Loaded interactively if not input.
%
%OUTPUT dataStruct : Same as input, with added field "frameAlignment".
%                    For each frame, it contains the subfields:
%           .centerOfMass: Center of mass of features in frame.
%           .eulerAnglesX: The Euler angles defining the rotation of this
%                          frame from its reference frame (which is the 4th
%                          entry in eulerAnglesX). Euler angles follow the
%                          x-convention. NaN means that no rotation was
%                          estimated for this frame, its coordinate system
%                          is obtained via the plane fit.
%           .coordSystem : Coordinate system of frame. If eulerAnglesX are
%                          NaN, it comes from the plane fit. If
%                          eulerAnglesX are not NaN, it is obtained by
%                          rotating the coordinate system of the reference
%                          frame.
%           .alignedCoord: Aligned coordinates in each frame, where both
%                          center of mass translation and rigid body
%                          rotation are compensated for.
%
%Khuloud Jaqaman, July 2007

%% preamble

%load dataStruct if not input
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end

%get number of frames in movie
numFrames = dataStruct.dataProperties.movieSize(end);

%put tracks into matrix format
[~,tracksIndx] = convStruct2MatNoMS(dataStruct.tracks);

%get total number of tracks
%numTracks = size(tracksIndx,1);

%assign 0 to outlier (unaligned & lagging) features in each frame
for iFrame = 1 : numFrames
    outIndx = [dataStruct.updatedClass(iFrame).unalignedIdx dataStruct.updatedClass(iFrame).laggingIdx];
    for iFeat = outIndx
        tracksIndx(tracksIndx(:,iFrame)==iFeat,iFrame) = 0;
    end
end

%reserve memory for euler angles and rotation frame of reference
eulerAnglesX = NaN(numFrames,4);

% old method
%
% %find which frames need their rotation to be estimated
% if isempty(dataStruct.planeFit)
%     frames2align = (1 : numFrames);
% else
%     frames2align = [];
%     for t = 1 : numFrames
%         if isempty(dataStruct.planeFit(t).planeVectors)
%             frames2align = [frames2align t];
%         end
%     end
% end
% framesOK = setxor((1:numFrames),frames2align);

% new method
if isempty(dataStruct.planeFit)
    frames2align = (1 : numFrames);
else
    % check for at least one plane-fit
    planeVectors = cat(1,dataStruct.planeFit.planeVectors);
    if isempty(planeVectors)
        frames2align = (1 : numFrames);
    else
        % get eigenValues
        eigenValues = catStruct(1,'dataStruct.planeFit.eigenValues');
        % find smallest eigen-value in each frame
        [~,minIdx] = min(eigenValues,[],2);
        % calc eigen ratios
        eigenRatios = NaN(numFrames,1);
        for iFrame = 1:numFrames
            eigenRatios(iFrame) = eigenValues(iFrame,minIdx(iFrame)) / mean(eigenValues(iFrame,setxor(1:3,minIdx(iFrame))));
        end
        % find frame(s) with smallest ratios
        [~,minFrame] = min(eigenRatios);
        frames2align = setxor(1:numFrames,minFrame);
    end
end
framesOK = setxor((1:numFrames),frames2align);

%if all frames are empty, make the first frame non-empty (it will be taken
%as the first reference frame)
if isempty(framesOK)
    frames2align = frames2align(2:end);
    framesOK = 1;
end

%find empty frames that are before the first non-empty frame
framesBefore1 = frames2align(frames2align < framesOK(1));
framesAfter1 = setxor(frames2align,framesBefore1); % the rest of the frames

%the reference frame of each empty frame before the first non-empty one is
%the frame after it
eulerAnglesX(framesBefore1,4) = framesBefore1 + 1;

%for the rest, their reference frame is the frame before each of them
eulerAnglesX(framesAfter1,4) = framesAfter1 - 1;

% %% center of mass shift
%
% %get the center of mass of each frame
% if ~isempty(dataStruct.planeFit)
%     centerOfMass = vertcat(dataStruct.planeFit.planeOrigin);
% else
%     centerOfMass = NaN(numFrames,3);
%     for iFrame = 1 : numFrames
%         centerOfMass(iFrame,:) = mean(dataStruct.initCoord(iFrame).allCoord(:,1:3));
%     end
% end
%
% %shift the track coordinates so that the origin is at the center of mass
% %in each frame
% for iFrame = 1 : numFrames
%     tracksInfo(:,8*(iFrame-1)+1:8*(iFrame-1)+3) = ...
%         tracksInfo(:,8*(iFrame-1)+1:8*(iFrame-1)+3) - ...
%         repmat(centerOfMass(iFrame,:),numTracks,1);
% end

%get the center of mass of each frame
if ~isempty(dataStruct.planeFit)
    centerOfMass = vertcat(dataStruct.planeFit.planeOrigin);
else
    centerOfMass = NaN(numFrames,3);
    for iFrame = 1 : numFrames
        centerOfMass(iFrame,:) = mean(dataStruct.initCoord(iFrame).allCoord(:,1:3));
    end
end

% get com corrected coords
allCoords(1:numFrames) = struct('coords',[]);
for iFrame = 1 : numFrames
    
    %get feature coordinates in this frame
    coordTmp = dataStruct.initCoord(iFrame).allCoord;
    numFeatures = size(coordTmp,1);
    
    %subtract center of mass from coordinates
    coordTmp(:,1:3) =  coordTmp(:,1:3) - repmat(centerOfMass(iFrame,:),numFeatures,1);
    
    % add to struct
    allCoords(iFrame).coords = coordTmp;
end

%% estimation of rotation

x0 = [0 0 0]; %initial guess
lb = [0 0 0]; %lower bound
ub = pi*[2 1 2]; %upper bound
options = optimset('Jacobian','on','Display','off','TolX', 1e-10, ...
    'Tolfun', 1e-10,'MaxFunEvals', 1e6, ...
    'MaxIter', 1e6); %minimization option - use analytical Jacobian

%reserve memory
coordSystem = repmat(eye(3),[1 1 numFrames]);

% %get the coordinate systems of non-empty frames from the plane fit
% for iFrame = framesOK
%     if isempty(dataStruct.planeFit) || isempty(dataStruct.planeFit(iFrame).planeVectors)
%         coordSystem(:,:,iFrame) = eye(3);
%     else
%         coordSystem(:,:,iFrame) = dataStruct.planeFit(iFrame).planeVectors;
%     end
% end

%record current coordinate systems
framesWPlane = [];
for iFrame = 1:numFrames
    if isempty(dataStruct.planeFit) || isempty(dataStruct.planeFit(iFrame).planeVectors)
        coordSystem(:,:,iFrame) = eye(3);
    else
        coordSystem(:,:,iFrame) = dataStruct.planeFit(iFrame).planeVectors;
        framesWPlane = [framesWPlane iFrame];
    end
end

framesNoPlane = setxor(1:numFrames,framesWPlane);

if ~isempty(framesWPlane)
    framesBefore = framesNoPlane(framesNoPlane < framesWPlane(1));
    framesAfter = framesNoPlane(framesNoPlane > framesWPlane(end));
    for iFrame = framesBefore(end:-1:1)
        coordSystem(:,:,iFrame) = coordSystem(:,:,iFrame+1);
    end
    for iFrame = framesAfter
        coordSystem(:,:,iFrame) = coordSystem(:,:,iFrame-1);
    end
end

%calculate the coordinate systems of empty frames that are before the first
%non-empty frame
%go from last to first and then frames after the last non empty frame
if ~isempty(frames2align)
    for iFrame = framesBefore1(end:-1:1)
        
        %determine reference frame
        jFrame = eulerAnglesX(iFrame,4);
        
        % new, EHarry, Nov 2012, get coords from coords and coordSystems
        [coord1, coord2] = planeCoords(allCoords,coordSystem,iFrame,jFrame,tracksIndx);
        
        % new, EHarry, Nov 2012, remove coordinates with unusually large or
        % small jumps using least median squares. Remove coords if either x
        % y or z jump is an outlier
        sep = coord1 - coord2;
        sep = sep';
        [~,~,~,outliers] = robustMean(sep(:),[],5);
        outliers = ceil(outliers/3);
        outliers = unique(outliers);
        coord1(outliers,:) = [];
        coord2(outliers,:) = [];
        
        % abort this frame if less than 3 coords remain
        if size(coord1,1) < 3
            continue
        end
        
        %estimate rotation angles (Euler angles following the x-convention)
        rotationAngles = lsqnonlin(@calcRotateResiduals,x0,lb,ub,options,coord1,coord2);
        
        %save rotation angles for output
        eulerAnglesX(iFrame,1:3) = rotationAngles;
        phi = eulerAnglesX(iFrame,1);
        theta = eulerAnglesX(iFrame,2);
        psi = eulerAnglesX(iFrame,3);
        
        %calculate the rotation matrix
        rotationMat1 = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
        rotationMat2 = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
        rotationMat3 = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
        rotationMat = rotationMat1 * rotationMat2 * rotationMat3;
        
        %rotate the coordinate system of this frame
        coordSystem(:,:,iFrame) = rotationMat \ coordSystem(:,:,iFrame);
        
    end
    for iFrame = framesAfter1
        
        %determine reference frame
        jFrame = eulerAnglesX(iFrame,4);
        
        % new, EHarry, Nov 2012, get coords from planeFit
        [coord1, coord2] = planeCoords(allCoords,coordSystem,iFrame,jFrame,tracksIndx);
        
        % new, EHarry, Nov 2012, remove coordinates with unusually large or
        % small jumps using least median squares. Remove coords if either x
        % y or z jump is an outlier
        sep = coord1 - coord2;
        sep = sep';
        [~,~,~,outliers] = robustMean(sep(:),[],5);
        outliers = ceil(outliers/3);
        outliers = unique(outliers);
        coord1(outliers,:) = [];
        coord2(outliers,:) = [];
        
        % abort this frame if less than 3 coords remain
        if size(coord1,1) < 3
            continue
        end
        
        %estimate rotation angles (Euler angles following the x-convention)
        rotationAngles = lsqnonlin(@calcRotateResiduals,x0,lb,ub,options,coord1,coord2);
        
        %save rotation angles for output
        eulerAnglesX(iFrame,1:3) = rotationAngles;
        phi = eulerAnglesX(iFrame,1);
        theta = eulerAnglesX(iFrame,2);
        psi = eulerAnglesX(iFrame,3);
        
        %calculate the rotation matrix
        rotationMat1 = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
        rotationMat2 = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
        rotationMat3 = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
        rotationMat = rotationMat1 * rotationMat2 * rotationMat3;
        
        %rotate the coordinate system of this frame
        coordSystem(:,:,iFrame) = rotationMat \ coordSystem(:,:,iFrame);
        
    end
    
end


%% coordinates in rotated & translated coordinate system

%reserve memory
alignedCoord(1:numFrames,1) = struct('values',[]);


%go over all frames
for iFrame = 1 : numFrames
    
    %     %get feature coordinates in this frame
    %     coordTmp = dataStruct.initCoord(iFrame).allCoord;
    %     numFeatures = size(coordTmp,1);
    %
    %     %subtract center of mass from coordinates
    %     coordTmp(:,1:3) =  coordTmp(:,1:3) - repmat(centerOfMass(iFrame,:),numFeatures,1);
    
    %get rotation matrix
    rotationMat = inv(coordSystem(:,:,iFrame));
    
    %calculate matrix to propagate error from original to rotated
    %coordinates
    errorPropMat = rotationMat.^2;
    
    %calculate coordinates in rotated coordinate system
    alignedCoord(iFrame).values(:,1:3) = (rotationMat * allCoords(iFrame).coords(:,1:3)')';
    
    %propagate error
    alignedCoord(iFrame).values(:,4:6) = sqrt((errorPropMat*((allCoords(iFrame).coords(:,4:6)).^2)')');
    
end




%% output to dataStruct

%fill in field
frameAlignment(1:numFrames,1) = struct('centerOfMass',[],'eulerAnglesX',[],...
    'coordSystem',[],'alignedCoord',[]);
for iFrame = 1 : numFrames
    frameAlignment(iFrame).centerOfMass = centerOfMass(iFrame,:);
    frameAlignment(iFrame).eulerAnglesX = eulerAnglesX(iFrame,:);
    frameAlignment(iFrame).coordSystem = coordSystem(:,:,iFrame);
    frameAlignment(iFrame).alignedCoord = alignedCoord(iFrame).values;
end


%% dependencies
% frameAlignment depends on dataProperites, initCoords, planeFit, tracks,
% sisterList and updatedClass
dependencies = struct('dataProperties',[],'initCoord',[],'planeFit',[],'tracks',[],'sisterList',[],'updatedClass',[]);

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

dependencies.dataProperties = dataPropV;
dependencies.initCoord = initCoordV;
dependencies.planeFit = planeFitV;
dependencies.tracks = tracksV;
dependencies.sisterList = sisterListV;
dependencies.updatedClass = updatedClassV;

% save into the first frameAlignment
frameAlignment(1).dependencies = dependencies;

%write field into dataStruct
dataStruct.frameAlignment = frameAlignment;


%% subfunction that calculates rotation residuals
    function [F,J] = calcRotateResiduals(x0,coord1,coord2)
        
        %fetch the Euler angles
        phi = x0(1);
        theta = x0(2);
        psi = x0(3);
        
        %calculate the rotation matrix
        rotationMatPsi = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
        rotationMatTheta = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
        rotationMatPhi = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
        rotationMat = rotationMatPsi * rotationMatTheta * rotationMatPhi;
        
        %calculate the residuals and put them in a 1D vector
        F = coord2 - coord1 * rotationMat';
        F = F(:);
        
        %calculate the Jacobian matrix - initialize
        J = zeros(length(F),3);
        
        %calculate the derivative of the rotation matrix with respect to phi
        rotationMatPhiD = [-sin(phi) cos(phi) 0; -cos(phi) -sin(phi) 0; 0 0 0];
        rotationMat = rotationMatPsi * rotationMatTheta * rotationMatPhiD;
        
        %calculate first column of Jacobian matrix
        JacobianCol = -coord1 * rotationMat';
        J(:,1) = JacobianCol(:);
        
        %calculate the derivative of the rotation matrix with respect to theta
        rotationMatThetaD = [0 0 0; 0 -sin(theta) cos(theta); 0 -cos(theta) -sin(theta)];
        rotationMat = rotationMatPsi * rotationMatThetaD * rotationMatPhi;
        
        %calculate second column of Jacobian matrix
        JacobianCol = -coord1 * rotationMat';
        J(:,2) = JacobianCol(:);
        
        %calculate the derivative of the rotation matrix with respect to psi
        rotationMatPsiD = [-sin(psi) cos(psi) 0; -cos(psi) -sin(psi) 0; 0 0 0];
        rotationMat = rotationMatPsiD * rotationMatTheta * rotationMatPhi;
        
        %calculate third column of Jacobian matrix
        JacobianCol = -coord1 * rotationMat';
        J(:,3) = JacobianCol(:);
    end

%% function to get plane coords of two tracks
    function [coord1, coord2] = planeCoords(allCoords,coordSystem,iFrame,jFrame,tracksIndx)
        %get the indices of features in this frame
        featIndx1 = tracksIndx(:,iFrame);
        
        %get the indices of features in reference frame
        featIndx2 = tracksIndx(:,jFrame);
        
        %find those features that are linked between the two frames
        goodIndx = featIndx1 ~= 0 & featIndx2 ~= 0;
        featIndx2 = featIndx2(goodIndx);
        featIndx1 = featIndx1(goodIndx);
        
        % make coords
        %coord1 = planeFit(jFrame).rotatedCoord(featIndx2,1:3);
        coord1 = (coordSystem(:,:,jFrame) \ allCoords(jFrame).coords(featIndx2,1:3)')';
        %coord2 = planeFit(iFrame).rotatedCoord(featIndx1,1:3);
        coord2 = (coordSystem(:,:,iFrame) \ allCoords(iFrame).coords(featIndx1,1:3)')';
    end
end
