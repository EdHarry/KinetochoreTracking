function dataStruct = makiPoleTransformSisters( dataStruct )

% EHarry Jan 2012

%% MAIN

if ~isfield(dataStruct,'poles');
    dataStruct.poleTransformedSisterList = [];
    return
end

if ~isfield(dataStruct,'sisterList');
    dataStruct.poleTransformedSisterList = [];
    return
end

if ~isfield(dataStruct,'planeFit');
    dataStruct.poleTransformedSisterList = [];
    return
end

if ~isfield(dataStruct,'tracks');
    dataStruct.poleTransformedSisterList = [];
    return
end

if isempty(dataStruct.poles)
    dataStruct.poleTransformedSisterList = [];
    return
end

if isempty(dataStruct.sisterList)
    dataStruct.poleTransformedSisterList = [];
    return
end

if isempty(dataStruct.planeFit)
    dataStruct.poleTransformedSisterList = [];
    return
end

if isempty(dataStruct.tracks)
    dataStruct.poleTransformedSisterList = [];
    return
end


poles = dataStruct.poles;
sisterList = dataStruct.sisterList;
tracks = dataStruct.tracks;
planeFit = dataStruct.planeFit;

pole1Track= poles.pole1Track;
pole2Track= poles.pole2Track;


frames=dataStruct.dataProperties.movieSize(4);

z=NaN(2*length(sisterList),frames);
theta=NaN(2*length(sisterList),frames);
r=NaN(2*length(sisterList),frames);
polePoleDistance = NaN(2,frames);
%k1 = NaN(length(sisterList),frames);
%k2 = NaN(length(sisterList),frames);
k1 = NaN(length(sisterList),3,frames,2);
k2 = NaN(length(sisterList),3,frames,2);
zeroD = NaN(1,frames);
pole1CoordsRotated = NaN(frames,3,2); % rotated coords of poles
pole2CoordsRotated = NaN(frames,3,2);
pole1Coords = NaN(frames,3,2); % original coords of poles
pole2Coords = NaN(frames,3,2);

for i = 1:length(pole1Track)
    pole1Complete_temp(:,:,i) = getTrackCoords(pole1Track(i) , frames);
    %pole1Complete_temp(:,:,i) = [pole1Track(i).tracksCoordAmpCG(1:8:end)',pole1Track(i).tracksCoordAmpCG(2:8:end)',pole1Track(i).tracksCoordAmpCG(3:8:end)',pole1Track(i).tracksCoordAmpCG(5:8:end)',pole1Track(i).tracksCoordAmpCG(6:8:end)',pole1Track(i).tracksCoordAmpCG(7:8:end)']; %complete pole tracks
end

for i = 1:length(pole2Track)
    pole2Complete_temp(:,:,i) = getTrackCoords(pole2Track(i) , frames);
    %pole2Complete_temp(:,:,i) = [pole2Track(i).tracksCoordAmpCG(1:8:end)',pole2Track(i).tracksCoordAmpCG(2:8:end)',pole2Track(i).tracksCoordAmpCG(3:8:end)',pole2Track(i).tracksCoordAmpCG(5:8:end)',pole2Track(i).tracksCoordAmpCG(6:8:end)',pole2Track(i).tracksCoordAmpCG(7:8:end)'];
end

pole1Complete(:,1:3) = nanmean(pole1Complete_temp(:,1:3,:),3);
pole2Complete(:,1:3) = nanmean(pole2Complete_temp(:,1:3,:),3);
pole1Complete(:,4:6) = sqrt(nansum(pole1Complete_temp(:,4:6,:).^2,3));
pole2Complete(:,4:6) = sqrt(nansum(pole2Complete_temp(:,4:6,:).^2,3));

polePoleDistance(1,:) = sqrt(sum(((pole1Complete(:,1:3) - pole2Complete(:,1:3)).^2),2));
polePoleDistance(2,:) = sqrt(sum(pole1Complete(:,4:6).^2,2)+sum(pole2Complete(:,4:6).^2,2));

for i = 1:length(pole1Track)
    featIdx1 = getFeatIdx(pole1Track(i),frames);
    featIdx2 = getFeatIdx(pole2Track(i),frames);
    
    for iTime=1:frames
        idX1 = featIdx1(iTime);
        idX2 = featIdx2(iTime);
        if ~isnan(idX1)
            pole1CoordsRotated(iTime,:,1,i) = planeFit(iTime).rotatedCoord(idX1,1:3);
            pole1CoordsRotated(iTime,:,2,i) = planeFit(iTime).rotatedCoord(idX1,4:6); % rotated coords of poles
        end
        if ~isnan(idX2)
            pole2CoordsRotated(iTime,:,1,i) = planeFit(iTime).rotatedCoord(idX2,1:3);
            pole2CoordsRotated(iTime,:,2,i) = planeFit(iTime).rotatedCoord(idX2,4:6);
        end
    end
end

iTrack=1;

for iTime = 1:frames
    pole1 = pole1Complete(iTrack,1:3);
    pole2 = pole2Complete(iTrack,1:3);
    pole1Error = pole1Complete(iTrack,4:6);
    pole2Error = pole2Complete(iTrack,4:6);
    
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
    
    iCyl=1;
    for j=1:length(sisterList) % convert to cylindrical coords
        if ~isempty(sisterList(j).coords1)
            if ~isnan(sisterList(j).coords1(iTime,1)) % check sister coords exist for timepoint
                
                s1Track = sisterList(1).trackPairs(j,1);
                s2Track = sisterList(1).trackPairs(j,2);
                
                %%%%%%%% V6 - go back to pole +ve and -ve pole definitions
                %             %%% need to swap sisters as well - see line 349
                %             s1Track = sisterList(1).trackPairs(j,2);
                %             s2Track = sisterList(1).trackPairs(j,1);
                %             %%%
                %%%%%%%%
%                 
%                 s1Complete = [tracks(s1Track).tracksCoordAmpCG(1:8:end)',tracks(s1Track).tracksCoordAmpCG(2:8:end)',tracks(s1Track).tracksCoordAmpCG(3:8:end)',tracks(s1Track).tracksCoordAmpCG(5:8:end)',tracks(s1Track).tracksCoordAmpCG(6:8:end)',tracks(s1Track).tracksCoordAmpCG(7:8:end)'];% complete sister tracks
%                 s2Complete = [tracks(s2Track).tracksCoordAmpCG(1:8:end)',tracks(s2Track).tracksCoordAmpCG(2:8:end)',tracks(s2Track).tracksCoordAmpCG(3:8:end)',tracks(s2Track).tracksCoordAmpCG(5:8:end)',tracks(s2Track).tracksCoordAmpCG(6:8:end)',tracks(s2Track).tracksCoordAmpCG(7:8:end)'];
%                 
%                 s1StartTime = tracks(s1Track).seqOfEvents(1,1);
%                 s1EndTime = tracks(s1Track).seqOfEvents(2,1); % find sister start and end of tracks
%                 
%                 s2StartTime = tracks(s2Track).seqOfEvents(1,1);
%                 s2EndTime = tracks(s2Track).seqOfEvents(2,1);
%                 
%                 sStartTime = max(s1StartTime,s2StartTime); % Take sister - siter start time
%                 sEndTime = min(s1EndTime,s2EndTime); % Take sister - siter end time
%                 
%                 s1TimePosition = iTime - s1StartTime +1 ;  % find where to find current sister time point
%                 s2TimePosition = iTime - s2StartTime +1 ;
%                 
%                 if iTime > sEndTime || iTime < sStartTime || s1TimePosition < 1  % grab sister position, leave a NaN if sister track extends beond pole track
%                     s1=[NaN,NaN,NaN];
%                 else
%                     s1 = s1Complete(s1TimePosition,1:3);
%                     s1Error = s1Complete(s1TimePosition,4:6);
%                 end
%                 
%                 if iTime > sEndTime || iTime < sStartTime || s2TimePosition < 1   % grab sister position, leave a NaN if sister track extends beond pole track
%                     s2=[NaN,NaN,NaN];
%                 else
%                     s2 = s2Complete(s2TimePosition,1:3);
%                     s2Error = s2Complete(s2TimePosition,4:6);
%                 end
%                 
%                 if isnan(s1(1)) || isnan(s2(1))  % leave only tracks present on both sisters
%                     s1=[NaN,NaN,NaN];
%                     s2=[NaN,NaN,NaN];
%                 end
                
                %  k1(j,iTime) = norm(s1-pole1);
                %  k2(j,iTime) = norm(s2-pole2);
                
                s1_full=getTrackCoords(tracks(s1Track) , frames);
                s2_full=getTrackCoords(tracks(s2Track) , frames);
                
                s1 = s1_full(iTime,1:3);
                s1Error = s1_full(iTime,4:6);
                
                s2 = s2_full(iTime,1:3);
                s2Error = s2_full(iTime,4:6);
                
                k1(j,1,iTime,1) = s1(1)-pole1(1);
                k1(j,2,iTime,1) = s1(2)-pole1(2);
                k1(j,3,iTime,1) = s1(3)-pole1(3);  % vectors between sisters and poles
                
                k2(j,1,iTime,1) = s2(1)-pole2(1);
                k2(j,2,iTime,1) = s2(2)-pole2(2);
                k2(j,3,iTime,1) = s2(3)-pole2(3);
                
                k1(j,1,iTime,2) = sqrt(s1Error(1).^2+pole1Error(1).^2);
                k1(j,2,iTime,2) = sqrt(s1Error(2).^2+pole1Error(2).^2);
                k1(j,3,iTime,2) = sqrt(s1Error(3).^2+pole1Error(3).^2);  % errors between sisters and poles
                
                k2(j,1,iTime,2) = sqrt(s2Error(1).^2+pole2Error(1).^2);
                k2(j,2,iTime,2) = sqrt(s2Error(2).^2+pole2Error(2).^2);
                k2(j,3,iTime,2) = sqrt(s2Error(3).^2+pole2Error(3).^2);
                
                z(iCyl,iTime) = sqrt((norm(s1-pole1))^2 - ((norm(cross(s1-pole1,s1-pole2)))^2/(norm(pole2-pole1))^2));   % z
                z(iCyl+1,iTime) = sqrt((norm(s2-pole1))^2 - ((norm(cross(s2-pole1,s2-pole2)))^2/(norm(pole2-pole1))^2));
                
                %%%%%%%% new for V6 - measure axis distances away from middle
                %%%%%%%% of spindle axis rather than a pole
                %%%%%%%% NOTE: this now makes this a LEFT-HAND coord system
                z(iCyl,iTime) = polePoleDistance(1,iTime)/2 - z(iCyl,iTime);
                z(iCyl+1,iTime) = polePoleDistance(1,iTime)/2 - z(iCyl+1,iTime);
                %%%%%%%%
                %%%%%%%%
                
                r(iCyl,iTime) = ((norm(cross(s1-pole1,s1-pole2)))/(norm(pole2-pole1)));   % r
                r(iCyl+1,iTime) = ((norm(cross(s2-pole1,s2-pole2)))/(norm(pole2-pole1)));
                
                
                
                %projection of kin onto p1 plane
                
                s1_p = s1 - [a,b,c]*((dot([a,b,c,d],[s1,1]))/(a^2+b^2+c^2));
                s2_p = s2 - [a,b,c]*((dot([a,b,c,d],[s2,1]))/(a^2+b^2+c^2));
                
                
                
                theta(iCyl,iTime) = atan2(dot(s1_p-pole1,refY_p),dot(s1_p-pole1,ref_p-pole1));   % theta
                if theta(iCyl,iTime) < 0
                    theta(iCyl,iTime) = 2*pi() - abs(theta(iCyl,iTime));
                end
                
                theta(iCyl+1,iTime) = atan2(dot(s2_p-pole1,refY_p),dot(s2_p-pole1,ref_p-pole1));
                if theta(iCyl+1,iTime) < 0
                    theta(iCyl+1,iTime) = 2*pi() - abs(theta(iCyl+1,iTime));
                end
                
                %          if acos((dot(s1_p-pole1,refY_p-pole1))/(norm(s1_p-pole1)*norm(refY_p-pole1))) < pi()/2
                %             theta(iCyl,iTime) = acos((dot(s1_p-pole1,ref_p-pole1))/(norm(s1_p-pole1)*norm(ref_p-pole1))); %theta
                %        else
                %           theta(iCyl,iTime) = 2*pi()-acos((dot(s1_p-pole1,ref_p-pole1))/(norm(s1_p-pole1)*norm(ref_p-pole1)));
                %      end
                %
                %      if acos((dot(s2_p-pole1,refY_p-pole1))/(norm(s2_p-pole1)*norm(refY_p-pole1))) < pi()/2
                %          theta(iCyl+1,iTime) = acos((dot(s2_p-pole1,ref_p-pole1))/(norm(s2_p-pole1)*norm(ref_p-pole1)));
                %      else
                %          theta(iCyl+1,iTime) = 2*pi()-acos((dot(s2_p-pole1,ref_p-pole1))/(norm(s2_p-pole1)*norm(ref_p-pole1)));
                %      end
                
                
            end
        end
        iCyl = iCyl+2;
    end
    
    iTrack=iTrack+1;
    
end

%% place everything in struct for neatness



poleTransformedSisterList.theta = theta;
poleTransformedSisterList.r = r;
poleTransformedSisterList.z = z;

poleTransformedSisterList.pole1CoordsRotated = pole1CoordsRotated;
poleTransformedSisterList.pole2CoordsRotated = pole2CoordsRotated;

dataStruct.poleTransformedSisterList = poleTransformedSisterList;

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


    function angles = getAngleWithNormal(trackCoords1,trackCoords2,normals)
        % takes the coords of two tracks and gives the angles (radians) of the vector
        % between then and the normal vectors
        
        trackVectors = trackCoords1 - trackCoords2; % take the vecotrs between the tracks
        
        [~,normals] = normList(normals); % normalise the normals
        [~,normedVectors]=normList(trackVectors); % normalise the track vecotrs
        
        angles = acos(dot(normedVectors,normals,2)); % calcule the angles
        
        angles(angles>pi/2) = pi - angles(angles>pi/2); % correct for incorrect direction of vectors
    end

    function normals = getNormals(planeFit)
        % gets the normals to the planeFit as a nx3 matrix, returns NaNs
        % for times with no normals
        
        normals = NaN(length(planeFit),3);
        
        for ii = 1:length(planeFit)
            vectors = planeFit(ii).planeVectors;
            if ~isempty(vectors)
                normals(ii,:) = vectors(:,1)';
            end
        end
    end

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

