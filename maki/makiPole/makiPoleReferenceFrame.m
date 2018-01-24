function dataStruct = makiPoleReferenceFrame( dataStruct,useAlignedCoords4Gaps,transformPoles )
%MAKIPOLETRANSFORMSISTERS Rotates all spots into a reference frame created
%by the spindle poles
% useAlignedCoords4Gaps ({1},0) - fills in missing frames with no poles with the
% relevent aligned coord system if it exists,
% EHarry Jan 2012

%% MAIN

if nargin < 2 || isempty(useAlignedCoords4Gaps)
    useAlignedCoords4Gaps = 1;
end

if nargin < 3 || isempty(transformPoles)
    transformPoles = 0;
end

if ~isfield(dataStruct,'poles');
    dataStruct.poleReferenceFrame = [];
    return
end

if ~isfield(dataStruct,'initCoord');
    dataStruct.poleReferenceFrame = [];
    return
end

if ~isfield(dataStruct,'frameAlignment') && useAlignedCoords4Gaps;
    dataStruct.poleReferenceFrame = [];
    return
end

if isempty(dataStruct.poles)
    dataStruct.poleReferenceFrame = [];
    return
end

if isempty(dataStruct.initCoord)
    dataStruct.poleReferenceFrame = [];
    return
end

if isempty(dataStruct.frameAlignment) && useAlignedCoords4Gaps
    dataStruct.poleReferenceFrame = [];
    return
end

poles = dataStruct.poles;
initCoord = dataStruct.initCoord;

if useAlignedCoords4Gaps
    frameAlignment = dataStruct.frameAlignment;
end

noPoles = [];

pole1Track= poles.pole1Track;
pole2Track= poles.pole2Track;


frames=dataStruct.dataProperties.movieSize(4);

poleReferenceFrame(1:frames) = struct('poleCoords_cartisian',[],'poleCoords_cylindrical',[],'origin',[],'planeVectors',[]);

for i = 1:length(pole1Track)
    pole1Complete_temp(:,:,i) = getTrackCoords(pole1Track(i) , frames);
    featIdx1(:,i) = getFeatIdx(pole1Track(i),frames);
end

for i = 1:length(pole2Track)
    pole2Complete_temp(:,:,i) = getTrackCoords(pole2Track(i) , frames);
    featIdx2(:,i) = getFeatIdx(pole2Track(i),frames);
end

pole1Complete(:,1:3) = nanmean(pole1Complete_temp(:,1:3,:),3);
pole2Complete(:,1:3) = nanmean(pole2Complete_temp(:,1:3,:),3);
pole1Complete(:,4:6) = sqrt(nansum(pole1Complete_temp(:,4:6,:).^2,3));
pole2Complete(:,4:6) = sqrt(nansum(pole2Complete_temp(:,4:6,:).^2,3));

polePoleDistance(:,1) = sqrt(sum(((pole1Complete(:,1:3) - pole2Complete(:,1:3)).^2),2));
polePoleDistance(:,2) = sqrt(sum(pole1Complete(:,4:6).^2,2)+sum(pole2Complete(:,4:6).^2,2));


for iTime = 1:frames
    
    poleSpotIdx = [featIdx1(iTime,:),featIdx2(iTime,:)];
    
    pole1 = pole1Complete(iTime,1:3);
    pole2 = pole2Complete(iTime,1:3);
    %     pole1Error = pole1Complete(iTime,4:6);
    %     pole2Error = pole2Complete(iTime,4:6);
    
    centre = (pole1+pole2)./2; % centre point between the two poles
    normal = pole1 - pole2;
    normal = normal./norm(normal); % normal vector between the poles
    planeVectors = calcPlaneVectors(normal'); % planeVectors of the coord system based on the poles
    
    a = pole2(1) - pole1(1);
    b = pole2(2) - pole1(2);
    c = pole2(3) - pole1(3);   % eqns of p1 plane
    
    d = -1*pole1(1)*(pole2(1)-pole1(1)) - pole1(2)*(pole2(2)-pole1(2)) - pole1(3)*(pole2(3)-pole1(3));
    
    % reference point
    
    ref = pole1 + [0,0,1];
    
    
    % projection of ref onto p1 plane
    
    ref_p = ref - [a,b,c]*((dot([a,b,c,d],[ref,1]))/(a^2+b^2+c^2));
    
    planeUnit = [a,b,c]'/(norm([a b c])); % unit normal to plane
    
    q0 = cos(pi()/4); % quaternion rotation of reference point around normal to plane by pi/2
    q1 = sin(pi()/4)*planeUnit(1);
    q2 = sin(pi()/4)*planeUnit(2);
    q3 = sin(pi()/4)*planeUnit(3);
    q = [(q0^2+q1^2-q2^2-q3^2),2*(q1*q2-q0*q3),2*(q1*q3+q0*q2);2*(q2*q1+q0*q3),(q0^2-q1^2+q2^2-q3^2),2*(q2*q3-q0*q1);2*(q3*q1-q0*q2),2*(q3*q2+q0*q1),(q0^2-q1^2-q2^2+q3^2)];
    
    refY_p = q*(ref_p-pole1)';
    refY_p = refY_p';
    
    allCoord = initCoord(iTime).allCoord;
    poleCoords_cylindrical = NaN(size(allCoord,1),3); % poleCoord_cyl go [z,r,theta]
    poleCoords_cartisian = NaN(size(allCoord,1),6); % poleCoord_cart go [x,y,z]
    
    if ~transformPoles
        spotsToTransform = find(~ismember(1:size(allCoord,1),poleSpotIdx)); % don't transform the poles
    else
        spotsToTransform = 1:size(allCoord,1); % if we want to transform the poles then transform everything
    end
    
    for iSpot = spotsToTransform
        if (~isnan(pole1(1)) && ~isnan(pole2(1))) || ~useAlignedCoords4Gaps % if the poles exist or we don't care about gaps
            spot = allCoord(iSpot,1:3);
            spotError = allCoord(iSpot,4:6);
            
            r = ((norm(cross(spot-pole1,spot-pole2)))/(norm(pole2-pole1)));   % r
            z = sqrt((norm(spot-pole1)).^2 - r.^2);   % z
            z = polePoleDistance(iTime,1)/2 - z; % this now makes this a LEFT-HAND coord system
            
            spot_p = spot - [a,b,c]*((dot([a,b,c,d],[spot,1]))/(a^2+b^2+c^2));
            
            theta = atan2(dot(spot_p-pole1,refY_p),dot(spot_p-pole1,ref_p-pole1));   % theta
            if theta < 0
                theta = 2*pi + theta;
            end
            
            poleCoords_cylindrical(iSpot,:) = [z,r,theta];
            
            %             [z,y,x] = pol2cart(theta,r,z); % transform back to cartisian coord relativd now to the static polar axis (note the rotation of x and z and the conversion to a right-handed system)
            %
            %             poleCoords_cartisian(iSpot,1:3) = [x,y,z];
            %             poleCoords_cartisian(iSpot,4:6) = sqrt(pole1Error.^2 + pole2Error.^2 + spotError.^2); % propogate the error through the cartisian transform, this overestimates the true error, TO DO: make this better by properly propergating the error into the cylindrical coordinates and then back out again
            
            rotationMat = inv(planeVectors); %rotation matrix
            errorPropMat = rotationMat.^2; %error propagation matrix
            poleCoords_cartisian(iSpot,1:3) = (rotationMat*(spot - centre)')';
            poleCoords_cartisian(iSpot,4:6) = sqrt((errorPropMat*((spotError).^2)')');
            
            
        elseif useAlignedCoords4Gaps % fill in gaps
            
            poleCoords_cylindrical(iSpot,:) = [NaN NaN NaN]; % don't fill in any cylindrical coords for this frame
            %poleCoords_cartisian(iSpot,:) = frameAlignment(iTime).alignedCoord(iSpot,:);
            poleCoords_cartisian(iSpot,:) = [NaN NaN NaN NaN NaN NaN]; % also don't fill in any spot info yet
            planeVectors = frameAlignment(iTime).coordSystem;
            %centre = frameAlignment(iTime).centerOfMass;
            centre = [NaN NaN NaN]; % don't fill in a centre yet
            noPoles = [noPoles iTime]; % update the list of frames with no poles
        end
    end
    
    %     poleCoords_cartisian = poleCoords_cartisian(~isnan(poleCoords_cartisian(:,1)),:); % eliminate the empty poles
    %     poleCoords_cylindrical = poleCoords_cylindrical(~isnan(poleCoords_cylindrical(:,1)),:);
    
    poleReferenceFrame(iTime).poleCoords_cylindrical = poleCoords_cylindrical;
    poleReferenceFrame(iTime).poleCoords_cartisian = poleCoords_cartisian;
    poleReferenceFrame(iTime).planeVectors = planeVectors;
    poleReferenceFrame(iTime).origin = centre;
end

% replcase the centres of the frames with no poles with a traslated centre
% of a frame with poles based on the traslation of the inlier spot com
% between the frames

if ~isempty(noPoles)
    noPoles = unique(noPoles); % frames with no poles, will be unique and sorted
    goodFrames = setxor(1:frames,noPoles); % good frames, sorted
    % frames before the first good frame will be targeted to the frame
    % after and procesed in decending order, frames after the first good frame will be targeted to the frame before and prcessed in acceding order
    before = noPoles(noPoles<goodFrames(1));
    after = noPoles(noPoles>goodFrames(1));
    before = sort(before,'descend');
    after = sort(after,'ascend');
    
    % for the frames before
    for frame = before
        target = frame + 1; % target the frame after
        traslationVector = frameAlignment(frame).centerOfMass - frameAlignment(target).centerOfMass; % vector from the target to the frame
        poleReferenceFrame(frame).origin = poleReferenceFrame(target).origin + traslationVector; % translate the origin
    end
    
    % for the frames after
    for frame = after
        target = frame - 1; % target the frame before
        traslationVector = frameAlignment(frame).centerOfMass - frameAlignment(target).centerOfMass; % vector from the target to the frame
        poleReferenceFrame(frame).origin = poleReferenceFrame(target).origin + traslationVector; % translate the origin
    end
    
    % now correct the spots in these frames by shifting them relative to the new centre
    for frame = noPoles
        allCoord = initCoord(frame).allCoord;
        poleSpotIdx = [featIdx1(frame,:),featIdx2(frame,:)];
        if ~transformPoles
            spotsToTransform = find(~ismember(1:size(allCoord,1),poleSpotIdx)); % don't transform the poles
        else
            spotsToTransform = 1:size(allCoord,1); % if we want to transform the poles then transform everything
        end
        rotationMat = inv(poleReferenceFrame(frame).planeVectors); %rotation matrix
        errorPropMat = rotationMat.^2; %error propagation matrix
        centre = poleReferenceFrame(frame).origin;
        %         spotTranslationVector = poleReferenceFrame(frame).origin - frameAlignment(frame).centerOfMass;
        %spotTranslationVector = frameAlignment(frame).centerOfMass - poleReferenceFrame(frame).origin;
        %poleReferenceFrame(frame).poleCoords_cartisian(:,1:3) = poleReferenceFrame(frame).poleCoords_cartisian(:,1:3) + repmat(spotTranslationVector,size(poleReferenceFrame(frame).poleCoords_cartisian(:,1),1),1);
        for iSpot = spotsToTransform
            spot = allCoord(iSpot,1:3);
            spotError = allCoord(iSpot,4:6);
            
            poleReferenceFrame(frame).poleCoords_cartisian(iSpot,1:3) = (rotationMat*(spot - centre)')';
            poleReferenceFrame(frame).poleCoords_cartisian(iSpot,4:6) = sqrt((errorPropMat*((spotError).^2)')');
        end
    end
    
    
end


%% dependencies
% poleReferenceFrame depends on dataProperites, initCoords, planeFit, tracks,
% sisterList, updatedClass, frameAlignment and poles
% dependencies = struct('dataProperties',[],'initCoord',[],'planeFit',[],'tracks',[],'sisterList',[],'updatedClass',[],'frameAlignment',[],'poles',[]);

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

% save into the first poleReferenceFrame
poleReferenceFrame(1).dependencies = dependencies;




dataStruct.poleReferenceFrame = poleReferenceFrame;


%%  SUBFUNTIONS

    function trackCoords = getTrackCoords(track , noTimePoints)
        % takes a track struct from maki and returns an nx3 matirx of the
        % track coords (absolute coords)
        
        trackCoords = NaN(noTimePoints,6); % predefine the coords, NaNs will be left for missing time points
        
        startTime = track.seqOfEvents(1,1); % get track start and end times
        endTime = track.seqOfEvents(2,1);
        
        X = track.tracksCoordAmpCG(1:8:end); % get the coords
        Y = track.tracksCoordAmpCG(2:8:end);
        Z = track.tracksCoordAmpCG(3:8:end);
        
        dX = track.tracksCoordAmpCG(5:8:end); % get the coords
        dY = track.tracksCoordAmpCG(6:8:end);
        dZ = track.tracksCoordAmpCG(7:8:end);
        
        trackCoords(startTime:endTime,:) = [X;Y;Z;dX;dY;dZ]';
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

