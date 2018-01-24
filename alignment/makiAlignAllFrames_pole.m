function dataStruct = makiAlignAllFrames_pole(dataStruct)
% edit of makiAlignFrames to align all frames to the frame with the "best"
% plane-fit, this is for pole based cordinates
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


if isempty(dataStruct.updatedClass_pole)
    dataStruct.frameAlignment_pole = [];
    return
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
if isempty(dataStruct.poleReferenceFrame)
    frames2align = (1 : numFrames);
else
    % get frames without poles
    framesWithoutPoles = [];
    for t = 1 : numFrames
        if all(isnan(dataStruct.poleReferenceFrame(t).poleCoords_cylindrical(:,1)))
            framesWithoutPoles = [framesWithoutPoles t];
        end
    end
    % get normals in z
    normalValues = catStruct(1,'dataStruct.poleReferenceFrame.planeVectors(3,3)');
    % set normals in frames without poles to Inf
    normalValues(framesWithoutPoles) = Inf;
    % find the frame wih the smallest z normal
    [~,minFrame] = min(abs(normalValues));
    % calc eigen ratios
    frames2align = setxor(1:numFrames,minFrame);
end
% frames2align = [];
%     for t = 1 : numFrames
%         if all(isnan(dataStruct.poleReferenceFrame(t).poleCoords_cylindrical(:,1)))
%             frames2align = [frames2align t];
%         end
%     end
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
centerOfMass = vertcat(dataStruct.poleReferenceFrame.origin);

% get com corrected coords
allCoords(1:numFrames) = struct('coords',[]);
for iFrame = 1 : numFrames
    
    %get feature coordinates in this frame
    coordTmp = dataStruct.initCoord(iFrame).allCoord;
    %numFeatures = size(coordTmp,1);
    
    %subtract center of mass from coordinates
    %coordTmp(:,1:3) =  coordTmp(:,1:3) - repmat(centerOfMass(iFrame,:),numFeatures,1);
    
    % add to struct
    allCoords(iFrame).coords = coordTmp;
end

%% estimation of rotation and shift

x0 = [0; 0; 0; 0; 0; 0]; %initial guess
%lb = [-Inf; -Inf; -Inf; 0; 0; 0]; %lower bound
%ub = [Inf; Inf; Inf; pi*2; pi; pi*2]; %upper bound
lb = [];
ub = [];
options = optimset('Jacobian','on','Display','off','TolX', 1e-10, ...
    'Tolfun', 1e-10,'MaxFunEvals', 1e6, ...
    'MaxIter', 1e6); %minimization option - use analytical Jacobian

%reserve memory
coordSystem = catStruct(3,'dataStruct.poleReferenceFrame.planeVectors');


%calculate the coordinate systems of empty frames that are before the first
%non-empty frame
%go from last to first and then frames after the last non empty frame
if ~isempty(frames2align)
    for iFrame = [framesBefore1(end:-1:1) framesAfter1]
        
        %determine reference frame
        jFrame = eulerAnglesX(iFrame,4);
        
        % new, EHarry, Nov 2012, get coords from coords and coordSystems
        [coord1, coord2] = planeCoords(allCoords,centerOfMass,coordSystem,iFrame,jFrame,tracksIndx);
        
        % abort this frame if less than 6 coords remain
        if size(coord1,1) < 6
            continue
        end
        
        % initial rotation guess
        rotationMat = coordSystem(:,:,jFrame) / coordSystem(:,:,iFrame);  
        [psi,theta,phi] = eulerAnglesFromRotMat(rotationMat);
        x0(4) = phi;
        x0(5) = theta;
        x0(6) = psi;
        
        %estimate rotation angles (Euler angles following the x-convention)
        rotationAngles = lsqnonlin(@calcRotateResiduals,x0,lb,ub,options,coord1,coord2);
        
        % save shift
        shift = rotationAngles(1:3)';
        
        %save rotation angles for output
        eulerAnglesX(iFrame,1:3) = rotationAngles(4:6)';
        phi = eulerAnglesX(iFrame,1);
        theta = eulerAnglesX(iFrame,2);
        psi = eulerAnglesX(iFrame,3);
        
        %calculate the rotation matrix
        rotationMat1 = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
        %rotationMatTheta = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
        rotationMat2 = [cos(theta) 0 sin(theta); 0 1 0;-sin(theta) 0 cos(theta)];
        %rotationMatPhi = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
        rotationMat3 = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
        rotationMat = rotationMat1 * rotationMat2 * rotationMat3;
        
        % take the shift away from the com
        centerOfMass(iFrame,:) = centerOfMass(iFrame,:) - shift;
        
        %rotate the coordinate system of this frame
        %coordSystem(:,:,iFrame) = (rotationMat' \ coordSystem(:,:,iFrame)')';
        coordSystem(:,:,iFrame) = rotationMat \ coordSystem(:,:,jFrame);
        
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
    
    % take away com from coords
    coords = allCoords(iFrame).coords(:,1:3) - repmat(centerOfMass(iFrame,:),size(allCoords(iFrame).coords(:,1:3),1),1);
    
    %calculate coordinates in rotated coordinate system
    alignedCoord(iFrame).values(:,1:3) = (rotationMat * coords')';
    
    %propagate error
    alignedCoord(iFrame).values(:,4:6) = sqrt((errorPropMat*((allCoords(iFrame).coords(:,4:6)).^2)')');
    
end




%% output to dataStruct

%fill in field
frameAlignment_pole(1:numFrames,1) = struct('centerOfMass',[],'eulerAnglesX',[],...
    'coordSystem',[],'alignedCoord',[]);
for iFrame = 1 : numFrames
    frameAlignment_pole(iFrame).centerOfMass = centerOfMass(iFrame,:);
    frameAlignment_pole(iFrame).eulerAnglesX = eulerAnglesX(iFrame,:);
    frameAlignment_pole(iFrame).coordSystem = coordSystem(:,:,iFrame);
    frameAlignment_pole(iFrame).alignedCoord = alignedCoord(iFrame).values;
end


%% dependencies
% frameAlignment_pole depends on dataProperites, initCoords, planeFit, tracks,
% sisterList, updatedClass, frameAlignment, poles, poleReferenceFrame,
% planeFit_pole, tracks_pole, sisterList_pole and updatedClass_pole
% dependencies = struct('dataProperties',[],'initCoord',[],'planeFit',[],'tracks',[],'sisterList',[],'updatedClass',[],'frameAlignment',[],'poles',[],'poleReferenceFrame',[],'planeFit_pole',[],'tracks_pole',[],'sisterList_pole',[],'updatedClass_pole',[]);

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
tracks_poleName = dataStruct.tracks_poleName;
tracks_poleV = getVersion(tracks_poleName);
sisterList_poleName = dataStruct.sisterList_poleName;
sisterList_poleV = getVersion(sisterList_poleName);
updatedClass_poleName = dataStruct.updatedClass_poleName;
updatedClass_poleV = getVersion(updatedClass_poleName);

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
dependencies.tracks_pole = tracks_poleV;
dependencies.sisterList_pole = sisterList_poleV;
dependencies.updatedClass_pole = updatedClass_poleV;

% save into the first frameAlignment_pole
frameAlignment_pole(1).dependencies = dependencies;

%write field into dataStruct
dataStruct.frameAlignment_pole = frameAlignment_pole;


%% subfunction that calculates rotation residuals
    function [F,J] = calcRotateResiduals(x0,coord1,coord2)
        
        % x0 -> [zShift;yShift;zShift;phi;theta;psi]
        
        % fetch shift
        shift = x0(1:3)';
        
        %fetch the Euler angles
        phi = x0(4);
        theta = x0(5);
        psi = x0(6);
        
        %calculate the rotation matrix
        rotationMatPsi = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
        %rotationMatTheta = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
        rotationMatTheta = [cos(theta) 0 sin(theta); 0 1 0;-sin(theta) 0 cos(theta)];
        %rotationMatPhi = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
        rotationMatPhi = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
        rotationMat = rotationMatPsi * rotationMatTheta * rotationMatPhi;
        
        %add the shift
        coord1 = coord1 + repmat(shift,size(coord1,1),1);
        
        %calculate the residuals and put them in a 1D vector
        F = coord1 * rotationMat' - coord2;
        F = F(:);
        
        %calculate the Jacobian matrix - initialize
        J = zeros(length(F),6);
        
        %calculate the derivative of the rotation matrix with respect to x
        %shift
        dx = zeros(size(coord1,1),3);
        dx(:,1) = 1;
        JacobianCol = dx * rotationMat';
        J(:,1) = JacobianCol(:);
        
        %calculate the derivative of the rotation matrix with respect to y
        %shift
        dy = zeros(size(coord1,1),3);
        dy(:,2) = 1;
        JacobianCol = dy * rotationMat';
        J(:,2) = JacobianCol(:);
        
        %calculate the derivative of the rotation matrix with respect to z
        %shift
        dz = zeros(size(coord1,1),3);
        dz(:,3) = 1;
        JacobianCol = dz * rotationMat';
        J(:,3) = JacobianCol(:);
        
        %calculate the derivative of the rotation matrix with respect to phi
        %rotationMatPhiD = [-sin(phi) cos(phi) 0; -cos(phi) -sin(phi) 0; 0 0 0];
        rotationMatPhiD = [0 0 0; 0 -sin(phi) -cos(phi); 0 cos(phi) -sin(phi)];
        rotationMat = rotationMatPsi * rotationMatTheta * rotationMatPhiD;
        
        %calculate first column of Jacobian matrix
        JacobianCol = coord1 * rotationMat';
        J(:,4) = JacobianCol(:);
        
        %calculate the derivative of the rotation matrix with respect to theta
        %rotationMatThetaD = [0 0 0; 0 -sin(theta) cos(theta); 0 -cos(theta) -sin(theta)];
        rotationMatThetaD = [-sin(theta) 0 cos(theta); 0 0 0; -cos(theta) 0 -sin(theta)];
        rotationMat = rotationMatPsi * rotationMatThetaD * rotationMatPhi;
        
        %calculate second column of Jacobian matrix
        JacobianCol = coord1 * rotationMat';
        J(:,5) = JacobianCol(:);
        
        %calculate the derivative of the rotation matrix with respect to psi
        rotationMatPsiD = [-sin(psi) -cos(psi) 0; cos(psi) -sin(psi) 0; 0 0 0];
        rotationMat = rotationMatPsiD * rotationMatTheta * rotationMatPhi;
        
        %calculate third column of Jacobian matrix
        JacobianCol = coord1 * rotationMat';
        J(:,6) = JacobianCol(:);
    end

%% function to get plane coords of two tracks
    function [coord1, coord2] = planeCoords(allCoords,com,coordSystem,iFrame,jFrame,tracksIndx)
        %get the indices of features in this frame
        featIndx1 = tracksIndx(:,iFrame);
        
        %get the indices of features in reference frame
        featIndx2 = tracksIndx(:,jFrame);
        
        %find those features that are linked between the two frames
        goodIndx = featIndx1 ~= 0 & featIndx2 ~= 0;
        featIndx2 = featIndx2(goodIndx);
        featIndx1 = featIndx1(goodIndx);
        
        % shift coords by com
        coord1 = allCoords(iFrame).coords(featIndx1,1:3);
        coord1 = coord1 - repmat(com(iFrame,:),size(coord1,1),1);
        coord2 = allCoords(jFrame).coords(featIndx2,1:3);
        coord2 = coord2 - repmat(com(jFrame,:),size(coord2,1),1);
        
        % make rotated coords
        coord1R = (coordSystem(:,:,iFrame) \ coord1')';
        coord2R = (coordSystem(:,:,jFrame) \ coord2')';
        
        % new, EHarry, Nov 2012, remove coordinates with unusually large or
        % small jumps using least median squares. Remove coords if either x
        % y or z jump is an outlier
        sep = coord1R - coord2R;
        sep = sep';
        [~,~,~,outliers] = robustMean(sep(:));
        outliers = ceil(outliers/3);
        outliers = unique(outliers);
        coord1(outliers,:) = [];
        coord2(outliers,:) = [];
    end
end
