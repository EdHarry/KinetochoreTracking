function [F,J] = modelTrackKT_objectiveFn( x0, param)
%MODELTRACKKT_OBJECTIVEFN objective function to be optimised by MODELTRACKKT
%   EHarry Dec 2012
% x0 -> [shift; scale; angles; ...; bleachingCorrection]
bleachingCorrection = x0(end);
shift = [x0(1:7:end-1) x0(2:7:end-1) x0(3:7:end-1)];
scale = x0(4:7:end-1);
angles = [x0(5:7:end-1) x0(6:7:end-1) x0(7:7:end-1)];
clear x0

% make image parameter matrix
imageParam = zeros(2*param.numPairs,8);
for iPair = 1:param.numPairs
    baseCoords = param.baseCoord(:,:,iPair);
    baseCentre = mean(baseCoords,2);
    phi = angles(iPair,1);
    theta = angles(iPair,2);
    psi = angles(iPair,3);
    %calculate the rotation matrix
    rotationMatPsi = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
    rotationMatTheta = [cos(theta) 0 sin(theta); 0 1 0;-sin(theta) 0 cos(theta)];
    rotationMatPhi = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    rotationMat = rotationMatPsi * rotationMatTheta * rotationMatPhi;
    % calculate translated coords
    coords = (rotationMat * (scale(iPair)*eye(3)) * (baseCoords - baseCentre(:,[1 1]))) + baseCentre(:,[1 1]) + repmat(shift(iPair,:)',1,2);
    idx = (iPair-1)*2 + 1;
    imageParam(idx:idx+1,1:3) = coords';
    imageParam(idx:idx+1,4) = bleachingCorrection * param.amp(iPair,:)';
end
imageParam(:,5) = param.psfSigma(1);
imageParam(:,6) = param.psfSigma(1);
imageParam(:,7) = param.psfSigma(2);
imageParam(:,8) = param.bg;

% fit parameters
fitParameters = [{'X1'} {'X2'} {'X3'} {'B'}];

% fit image
[~,~,~,~,~,~,resAndGauss,jacobianXYZ] = GaussFitND_noFit(param.image,[],fitParameters,imageParam);

% set F
F = resAndGauss(:,1);

% make J
J = zeros(length(F),(7*param.numPairs)+1);
for iPair = 1:param.numPairs
    idx = (iPair-1)*7 + 1;
    idxXYZ = (iPair-1)*6 + 1;
    % jacobian rel to shift = sum jabobain x y and z
    J(:,idx:idx+2) = repmat(sum(jacobianXYZ(:,idxXYZ:idxXYZ+5),2),1,3);
    % jacobian rel to scale is (sum jabobain x y and z) .* rotmat * shifted coords
    rotCoord = rotationMat * (baseCoords - baseCentre(:,[1 1]));
    J(:,idx+3) = jacobianXYZ(:,idxXYZ).*rotCoord(1,1) + jacobianXYZ(:,idxXYZ+1).*rotCoord(2,1) + jacobianXYZ(:,idxXYZ+2).*rotCoord(3,1) + jacobianXYZ(:,idxXYZ+3).*rotCoord(1,2) + jacobianXYZ(:,idxXYZ+4).*rotCoord(2,2) + jacobianXYZ(:,idxXYZ+5).*rotCoord(3,2);
    % jacobian rel to angles, make dif rot mat and set for each angle
    % wrt phi
    rotationMatPhiD = [0 0 0; 0 -sin(phi) -cos(phi); 0 cos(phi) -sin(phi)];
    rotationMat = rotationMatPsi * rotationMatTheta * rotationMatPhiD;
    rotCoord = rotationMat * (baseCoords - baseCentre(:,[1 1]));
    J(:,idx+4) = jacobianXYZ(:,idxXYZ).*rotCoord(1,1) + jacobianXYZ(:,idxXYZ+1).*rotCoord(2,1) + jacobianXYZ(:,idxXYZ+2).*rotCoord(3,1) + jacobianXYZ(:,idxXYZ+3).*rotCoord(1,2) + jacobianXYZ(:,idxXYZ+4).*rotCoord(2,2) + jacobianXYZ(:,idxXYZ+5).*rotCoord(3,2);
    % wrt theta
    rotationMatThetaD = [-sin(theta) 0 cos(theta); 0 0 0; -cos(theta) 0 -sin(theta)];
    rotationMat = rotationMatPsi * rotationMatThetaD * rotationMatPhi;
    rotCoord = rotationMat * (baseCoords - baseCentre(:,[1 1]));
    J(:,idx+5) = jacobianXYZ(:,idxXYZ).*rotCoord(1,1) + jacobianXYZ(:,idxXYZ+1).*rotCoord(2,1) + jacobianXYZ(:,idxXYZ+2).*rotCoord(3,1) + jacobianXYZ(:,idxXYZ+3).*rotCoord(1,2) + jacobianXYZ(:,idxXYZ+4).*rotCoord(2,2) + jacobianXYZ(:,idxXYZ+5).*rotCoord(3,2);
    % wrt psi
    rotationMatPsiD = [-sin(psi) -cos(psi) 0; cos(psi) -sin(psi) 0; 0 0 0];
    rotationMat = rotationMatPsiD * rotationMatTheta * rotationMatPhi;
    rotCoord = rotationMat * (baseCoords - baseCentre(:,[1 1]));
    J(:,idx+6) = jacobianXYZ(:,idxXYZ).*rotCoord(1,1) + jacobianXYZ(:,idxXYZ+1).*rotCoord(2,1) + jacobianXYZ(:,idxXYZ+2).*rotCoord(3,1) + jacobianXYZ(:,idxXYZ+3).*rotCoord(1,2) + jacobianXYZ(:,idxXYZ+4).*rotCoord(2,2) + jacobianXYZ(:,idxXYZ+5).*rotCoord(3,2);
end
J(:,end) = resAndGauss(:,2) ./ bleachingCorrection;

end

