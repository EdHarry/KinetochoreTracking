function [eflag,poleCoord] = cylindericalTransform(planeFit,tracks,sisterList)
%poles find poles
%   find poles in centrin movies and convert sisters to spindle axis coords
% EHarry October 2011

eflag=0; % no poles found flag
poleCoord=0;

%%% pole cirteria
maxDistance=12; % max pole search distance (microns)
minDistance=5; % min pole search distacne (microns)
maxAngle=50*pi()/180; % max pole search angle (degrees in radians)
minFrames=10; % min no. of frames a pole track must exist
%%%

foundPositivePole=0; % found pole markers
foundNegativePole=0;

frames=length(planeFit);

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

%%%%%%%%%%%%% Old Stuff

%planeOrign = NaN(3,frames);
%planeOrignRaw = planeOrign;

%for i=1:frames
%    planeOrignRaw(:,i) = planeFit(i).planeOrigin;
%end

%     for iTime=1:frames-10 % loop until poles are found
%
%         unInd = planeFit(iTime).unalignedIdx; % unaligned spot index
%         unXYZ = planeFit(iTime).rotatedCoord(unInd,1:3); % XYZ coords of unaligned spots
%
%         unIndPositive = unInd(unXYZ(:,1)>0); % index of unaligned spots to the +ve side of the plate
%         unIndNegative = unInd(unXYZ(:,1)<0); % index of unaligned spots to the -ve side of the plate
%         unXYZPositive = unXYZ(unXYZ(:,1)>0,:); % XYZ coords of unaligned spots to the +ve side of the plate
%         unXYZNegative = unXYZ(unXYZ(:,1)<0,:); % XYZ coords of unaligned spots to the -ve side of the plate
%
%         unXYZPositive_5 = unXYZPositive(abs(unXYZPositive(:,1))>5,:); % filter out spots less than 5 micros away from the plate
%         unXYZNegative_5 = unXYZNegative(abs(unXYZNegative(:,1))>5,:);
%         unIndPositive_5 = unIndPositive(abs(unXYZPositive(:,1))>5);
%         unIndNegative_5 = unIndNegative(abs(unXYZNegative(:,1))>5);
%
%         unXYZPositive = unXYZPositive_5(abs(unXYZPositive_5(:,1))<8,:); % filter out spots more than 8 micros away from the plate
%         unXYZNegative = unXYZNegative_5(abs(unXYZNegative_5(:,1))<8,:);
%         unIndPositive = unIndPositive_5(abs(unXYZPositive_5(:,1))<8);
%         unIndNegative = unIndNegative_5(abs(unXYZNegative_5(:,1))<8);
%
%
%
%
%
%            for i=1:size(unXYZPositive,1)   % loop over all posible pairs of pole candidates
%                 for j=1:size(unXYZNegative,1)
%                      foundPositivePole=0;
%                      foundNegativePole=0;
%
%                    % angle = acos((dot(unXYZPositive(i,:)-unXYZNegative(j,:),[1,0,0]))/(norm(unXYZPositive(i,:)-unXYZNegative(j,:))));  % calculate angle of candidates with rotated plane fit
%                    angle1 = acos(dot([1,0,0],abs(unXYZPositive(i,:)))/norm(abs(unXYZPositive(i,:)))); %% changed to calculate angle of each candidate with plate COM (0,0,0)
%                    angle2 = acos(dot([1,0,0],abs(unXYZNegative(j,:)))/norm(abs(unXYZNegative(j,:))));
%
%                    if angle1 < 15*pi()/180 && angle2 < 15*pi()/180 % check that both angles are no more than 15 degrees
%
%                            candIndPositive = unIndPositive(i); % index of possible candidates
%                            candIndNegative = unIndNegative(j);
%
%                            for iTrack = 1:length(tracks) % find the tracks of the candidates
%                                 if tracks(iTrack).seqOfEvents(1,1) == iTime % check that track begins at the current frame
%
%                                     if tracks(iTrack).tracksFeatIndxCG(1) == candIndPositive % find positive candidate track
%                                         if tracks(iTrack).seqOfEvents(2,1) - tracks(iTrack).seqOfEvents(1,1) + 1 > 10 % check that track is more than 10 frames long
%                                             pole1Track = iTrack; % found pole1
%                                             foundPositivePole=1;
%
%                                         end
%
%                                     elseif tracks(iTrack).tracksFeatIndxCG(1) == candIndNegative % find negative candidate track
%                                         if tracks(iTrack).seqOfEvents(2,1) - tracks(iTrack).seqOfEvents(1,1) + 1 > 10 % check that track is more than 10 frames long
%                                             pole2Track = iTrack; % found pole2
%                                             foundNegativePole=1;
%
%                                         end
%
%
%                                     end
%                                 end
%                            end
%                     end
%
%                     if foundPositivePole==1 && foundNegativePole==1 % check whether both poles have been found
%                         trackStart=iTime;
%                         polesFound=1;
%                         break
%                     end
%
%                 end
%
%                 if polesFound==1
%                     break
%                 end
%
%            end
%
%
%
%
%
%
%         if polesFound==1
%            break
%         end
%
%
%     end

%%%%%%%%%%%%%%%%%%

%%% New more robust method, look for each pole on its own, and check citeria each frame of the track

%%% pole on positive side (pole1)
for iTime=1:frames-minFrames % loop until poles are found
    
    unInd = planeFit(iTime).unalignedIdx; % unaligned spot index
    unXYZ = planeFit(iTime).rotatedCoord(unInd,1:3); % XYZ coords of unaligned spots
    
    unIndPositive = unInd(unXYZ(:,1)>0); % index of unaligned spots to the +ve side of the plate
    unXYZPositive = unXYZ(unXYZ(:,1)>0,:); % XYZ coords of unaligned spots to the +ve side of the plate
    
    unXYZPositive_min = unXYZPositive(abs(unXYZPositive(:,1))>minDistance,:); % filter out spots less than minDistance microns away from the plate
    unIndPositive_min = unIndPositive(abs(unXYZPositive(:,1))>minDistance);
    
    unXYZPositive = unXYZPositive_min(abs(unXYZPositive_min(:,1))<maxDistance,:); % filter out spots more than maxDistance microns away from the plate
    unIndPositive = unIndPositive_min(abs(unXYZPositive_min(:,1))<maxDistance);
    
    for i=1:size(unXYZPositive,1)   % loop over all posible pairs of pole candidates
        
        foundPositivePole=0;
        
        
        angle1 = acos(dot([1,0,0],abs(unXYZPositive(i,:)))/norm(abs(unXYZPositive(i,:)))); %% changed to calculate angle of each candidate with plate COM (0,0,0)
        
        
        if angle1 < maxAngle  % check that angle is no more than maxAngle
            
            candIndPositive = unIndPositive(i); % index of possible candidates
            
            for iTrack = 1:length(tracks) % find the tracks of the candidates
                if tracks(iTrack).seqOfEvents(1,1) == iTime % check that track begins at the current frame
                    
                    if tracks(iTrack).tracksFeatIndxCG(1) == candIndPositive % find positive candidate track
                        if tracks(iTrack).seqOfEvents(2,1) - tracks(iTrack).seqOfEvents(1,1) + 1 >= minFrames % check that track is more than minFrames long
                            
                            goodTrack=1; % assume track is good unless it failes a condtion test in any of the next frames
                            trackIdx=2; % track indx, starting at second track point
                            for iTime2 = iTime+1:tracks(iTrack).seqOfEvents(2,1) % check rest of the frame in the movie that the poleTrack stays in the right place
                                
                                featIdx = tracks(iTrack).tracksFeatIndxCG(trackIdx); % spot index in next frame to test
                                testUnIdx = planeFit(iTime2).unalignedIdx;
                                if featIdx~=0 % if a track spot exists in the next frame
                                    if sum(featIdx==testUnIdx)==1 % check in next frame that the spot is still unaligned
                                        XYZ = planeFit(iTime2).rotatedCoord(featIdx,1:3); % rotated coords of next spot
                                        if abs(XYZ(1)) > minDistance && abs(XYZ(1)) < maxDistance % check spot is still between minDistance and maxDistance microns
                                            angle = acos(dot([1,0,0],abs(XYZ))/norm(abs(XYZ)));
                                            if angle < maxAngle % check angle is still less than maxAngle degrees
                                            else
                                                goodTrack=0; % poleTrack is not good, indicate as such and break loop
                                                break
                                            end
                                        else
                                            goodTrack=0; % poleTrack is not good, indicate as such and break loop
                                            break
                                        end
                                    else
                                        goodTrack=0; % poleTrack is not good, indicate as such and break loop
                                        break
                                    end
                                end
                                trackIdx=trackIdx+1; % inc trackInx
                            end
                            
                            if goodTrack==1 % if track is still good
                                pole1Track = iTrack;
                                foundPositivePole=1; % indicate the pole has been found
                            end
                            
                        end
                        
                        
                        
                    end
                end
                
                if foundPositivePole==1 % check whether pole has been found and exit loop early if it has
                    break
                end
                
            end
        end
        
        if foundPositivePole==1  % check whether pole has been found and exit loop early if it has
            break
        end
        
        
    end
    
    
    
    
    
    
    if foundPositivePole==1
        trackStartPositive=iTime;
        break
    end
    
    
end

%%%

%%% pole on negative side (pole2)
for iTime=1:frames-10 % loop until poles are found
    
    unInd = planeFit(iTime).unalignedIdx; % unaligned spot index
    unXYZ = planeFit(iTime).rotatedCoord(unInd,1:3); % XYZ coords of unaligned spots
    
    unIndNegative = unInd(unXYZ(:,1)<0); % index of unaligned spots to the -ve side of the plate
    unXYZNegative = unXYZ(unXYZ(:,1)<0,:); % XYZ coords of unaligned spots to the -ve side of the plate
    
    unXYZNegative_min = unXYZNegative(abs(unXYZNegative(:,1))>minDistance,:); % filter out spots less than minDistance microns away from the plate
    unIndNegative_min = unIndNegative(abs(unXYZNegative(:,1))>minDistance);
    
    unXYZNegative = unXYZNegative_min(abs(unXYZNegative_min(:,1))<maxDistance,:); % filter out spots more than maxDistance microns away from the plate
    unIndNegative = unIndNegative_min(abs(unXYZNegative_min(:,1))<maxDistance);
    
    for i=1:size(unXYZNegative,1)   % loop over all posible pairs of pole candidates
        
        foundNegativePole=0;
        
        
        angle2 = acos(dot([1,0,0],abs(unXYZNegative(i,:)))/norm(abs(unXYZNegative(i,:)))); %% changed to calculate angle of each candidate with plate COM (0,0,0)
        
        
        if angle2 < maxAngle  % check that angle is no more than maxAngle degrees
            
            candIndNegative = unIndNegative(i); % index of possible candidates
            
            for iTrack = 1:length(tracks) % find the tracks of the candidates
                if tracks(iTrack).seqOfEvents(1,1) == iTime % check that track begins at the current frame
                    
                    if tracks(iTrack).tracksFeatIndxCG(1) == candIndNegative % find positive candidate track
                        if tracks(iTrack).seqOfEvents(2,1) - tracks(iTrack).seqOfEvents(1,1) + 1 >= minFrames % check that track is more than minFrames long
                            
                            goodTrack=1; % assume track is good unless it failes a condtion test in any of the next frames
                            trackIdx=2; % track indx, starting at second track point
                            for iTime2 = iTime+1:tracks(iTrack).seqOfEvents(2,1) % check rest of the frame in the movie that the poleTrack stays in the right place
                                
                                featIdx = tracks(iTrack).tracksFeatIndxCG(trackIdx); % spot index in next frame to test
                                testUnIdx = planeFit(iTime2).unalignedIdx;
                                if featIdx~=0 % if a track spot exists in the next frame
                                    if sum(featIdx==testUnIdx)==1 % check in next frame that the spot is still unaligned
                                        XYZ = planeFit(iTime2).rotatedCoord(featIdx,1:3); % rotated coords of next spot
                                        if abs(XYZ(1)) > minDistance && abs(XYZ(1)) < maxDistance % check spot is still between minDistance and maxDistance microns
                                            angle = acos(dot([1,0,0],abs(XYZ))/norm(abs(XYZ)));
                                            if angle < maxAngle % check angle is still less than maxAngle
                                            else
                                                goodTrack=0; % poleTrack is not good, indicate as such and break loop
                                                break
                                            end
                                        else
                                            goodTrack=0; % poleTrack is not good, indicate as such and break loop
                                            break
                                        end
                                    else
                                        goodTrack=0; % poleTrack is not good, indicate as such and break loop
                                        break
                                    end
                                end
                                trackIdx=trackIdx+1; % inc trackInx
                            end
                            
                            if goodTrack==1 % if track is still good
                                pole2Track = iTrack;
                                foundNegativePole=1; % indicate the pole has been found
                            end
                            
                        end
                        
                        
                        
                    end
                end
                
                if foundNegativePole==1 % check whether pole has been found and exit loop early if it has
                    break
                end
                
            end
        end
        
        if foundNegativePole==1  % check whether pole has been found and exit loop early if it has
            break
        end
        
        
    end
    
    
    
    
    
    
    if foundNegativePole==1
        trackStartNegative=iTime;
        break
    end
    
    
end

%%%

if foundPositivePole==0 || foundNegativePole==0
    eflag=1;
    return
end



pole1Complete = [tracks(pole1Track).tracksCoordAmpCG(1:8:end)',tracks(pole1Track).tracksCoordAmpCG(2:8:end)',tracks(pole1Track).tracksCoordAmpCG(3:8:end)',tracks(pole1Track).tracksCoordAmpCG(5:8:end)',tracks(pole1Track).tracksCoordAmpCG(6:8:end)',tracks(pole1Track).tracksCoordAmpCG(7:8:end)']; %complete pole tracks
pole2Complete = [tracks(pole2Track).tracksCoordAmpCG(1:8:end)',tracks(pole2Track).tracksCoordAmpCG(2:8:end)',tracks(pole2Track).tracksCoordAmpCG(3:8:end)',tracks(pole2Track).tracksCoordAmpCG(5:8:end)',tracks(pole2Track).tracksCoordAmpCG(6:8:end)',tracks(pole2Track).tracksCoordAmpCG(7:8:end)'];


%%%%%%%%%%%%%% V6 - back to old pole definitions
% %%%% swap poles to make data easier to compare to old coords
%
% temp = pole1Complete;
% pole1Complete = pole2Complete;
% pole2Complete = temp;
%
% temp = pole1Track;
% pole1Track = pole2Track;
% pole2Track = temp;
%
% temp = trackStartPositive;
% trackStartPositive = trackStartNegative;
% trackStartNegative = temp;
%
% %%%%%%
%%%%%%%%%%%%%%


%%convert to cylindrical coords

trackStart = max(trackStartPositive,trackStartNegative); % latest start frame

pole1Complete = pole1Complete(trackStart-trackStartPositive+1:end,:); % cut out early track coords with only one pole
pole2Complete = pole2Complete(trackStart-trackStartNegative+1:end,:);

trackEnd = min(size(pole1Complete,1),size(pole2Complete,1)) + trackStart -1; % earliest end frame

pole1Complete = pole1Complete(1:min(size(pole1Complete,1),size(pole2Complete,1)),:); % cut out late track coords with only one pole
pole2Complete = pole2Complete(1:min(size(pole1Complete,1),size(pole2Complete,1)),:);

polePoleD = sqrt(sum(((pole1Complete(:,1:3) - pole2Complete(:,1:3)).^2),2));
polePoleError = sqrt(sum(pole1Complete(:,4:6).^2,2)+sum(pole2Complete(:,4:6).^2,2));

polePoleDistance(1,trackStart:trackEnd) = polePoleD; % pole pole distance
polePoleDistance(2,trackStart:trackEnd) = polePoleError;

pole1Coords(trackStart:trackEnd,:,1) = pole1Complete(:,1:3); % coords of poles
pole2Coords(trackStart:trackEnd,:,1) = pole2Complete(:,1:3);
pole1Coords(trackStart:trackEnd,:,2) = pole1Complete(:,4:6); % coords of poles
pole2Coords(trackStart:trackEnd,:,2) = pole2Complete(:,4:6);


%track spot indexs
featInx1 = tracks(pole1Track).tracksFeatIndxCG;
featInx2 = tracks(pole2Track).tracksFeatIndxCG;

featInx1=featInx1(trackStart-trackStartPositive+1:end);
featInx2=featInx2(trackStart-trackStartNegative+1:end);

featInx1=featInx1(1:min(length(featInx1),length(featInx2)));
featInx2=featInx2(1:min(length(featInx1),length(featInx2)));

idX=1; % index counter
for iTime=trackStart:trackEnd
    idX1 = featInx1(idX);
    idX2 = featInx2(idX);
    if idX1~=0
        pole1CoordsRotated(iTime,:,1) = planeFit(iTime).rotatedCoord(idX1,1:3);
        pole1CoordsRotated(iTime,:,2) = planeFit(iTime).rotatedCoord(idX1,4:6); % rotated coords of poles
    end
    if idX2~=0
        pole2CoordsRotated(iTime,:,1) = planeFit(iTime).rotatedCoord(idX2,1:3);
        pole2CoordsRotated(iTime,:,2) = planeFit(iTime).rotatedCoord(idX2,4:6);
    end
    idX=idX+1;
end


iTrack=1;


%%% for planefit distance away from pole
for iTime = trackStart:trackEnd
    
    o4 = pole1Complete(iTrack,1:3);
    o5 = pole2Complete(iTrack,1:3);
    o=planeFit(iTime).planeOrigin;
    plane=planeFit(1).plane;
    o2(1)=o(1);
    o2(2)=o(2)+1;
    o2(3)=(plane(4)-plane(1)*o2(1)-plane(2)*o2(2))/plane(3);
    o3(1)=o(1)+1;
    o3(2)=o(2);
    o3(3)=(plane(4)-plane(1)*o3(1)-plane(2)*o3(2))/plane(3);
    
    t=-det([1,1,1,1 ; o(1),o2(1),o3(1),o4(1) ; o(2),o2(2),o3(2),o4(2) ; o(3),o2(3),o3(3),o4(3)])/det([1,1,1,0 ; o(1),o2(1),o3(1),o5(1)-o4(1) ; o(2),o2(2),o3(2),o5(2)-o4(2) ; o(3),o2(3),o3(3),o5(3)-o4(3)]);
    
    o6=[o4(1)+t*(o5(1)-o4(1)),o4(2)+t*(o5(2)-o4(2)),o4(3)+t*(o5(3)-o4(3))];
    
    zeroD(iTime) = sqrt(sum(((o4 - o6).^2),2));
    
    iTrack = iTrack+1;
end

iTrack=1;

for iTime = trackStart:trackEnd
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
                
                s1Complete = [tracks(s1Track).tracksCoordAmpCG(1:8:end)',tracks(s1Track).tracksCoordAmpCG(2:8:end)',tracks(s1Track).tracksCoordAmpCG(3:8:end)',tracks(s1Track).tracksCoordAmpCG(5:8:end)',tracks(s1Track).tracksCoordAmpCG(6:8:end)',tracks(s1Track).tracksCoordAmpCG(7:8:end)'];% complete sister tracks
                s2Complete = [tracks(s2Track).tracksCoordAmpCG(1:8:end)',tracks(s2Track).tracksCoordAmpCG(2:8:end)',tracks(s2Track).tracksCoordAmpCG(3:8:end)',tracks(s2Track).tracksCoordAmpCG(5:8:end)',tracks(s2Track).tracksCoordAmpCG(6:8:end)',tracks(s2Track).tracksCoordAmpCG(7:8:end)'];
                
                s1StartTime = tracks(s1Track).seqOfEvents(1,1);
                s1EndTime = tracks(s1Track).seqOfEvents(2,1); % find sister start and end of tracks
                
                s2StartTime = tracks(s2Track).seqOfEvents(1,1);
                s2EndTime = tracks(s2Track).seqOfEvents(2,1);
                
                sStartTime = max(s1StartTime,s2StartTime); % Take sister - siter start time
                sEndTime = min(s1EndTime,s2EndTime); % Take sister - siter end time
                
                s1TimePosition = iTime - s1StartTime +1 ;  % find where to find current sister time point
                s2TimePosition = iTime - s2StartTime +1 ;
                
                if iTime > sEndTime || iTime < sStartTime || s1TimePosition < 1  % grab sister position, leave a NaN if sister track extends beond pole track
                    s1=[NaN,NaN,NaN];
                else
                    s1 = s1Complete(s1TimePosition,1:3);
                    s1Error = s1Complete(s1TimePosition,4:6);
                end
                
                if iTime > sEndTime || iTime < sStartTime || s2TimePosition < 1   % grab sister position, leave a NaN if sister track extends beond pole track
                    s2=[NaN,NaN,NaN];
                else
                    s2 = s2Complete(s2TimePosition,1:3);
                    s2Error = s2Complete(s2TimePosition,4:6);
                end
                
                if isnan(s1(1)) || isnan(s2(1))  % leave only tracks present on both sisters
                    s1=[NaN,NaN,NaN];
                    s2=[NaN,NaN,NaN];
                end
                
                %  k1(j,iTime) = norm(s1-pole1);
                %  k2(j,iTime) = norm(s2-pole2);
                
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

clear poleCoord % clear poleCoord=0 to turn into struct

poleCoord.theta = theta;
poleCoord.r = r;
poleCoord.z = z;
poleCoord.k1 = k1;
poleCoord.k2 = k2;
poleCoord.polePoleDistance = polePoleDistance;
poleCoord.zeroD = zeroD;
poleCoord.pole1CoordsRotated = pole1CoordsRotated;
poleCoord.pole2CoordsRotated = pole2CoordsRotated;
poleCoord.pole1Coords = pole1Coords;
poleCoord.pole2Coords = pole2Coords;

end

