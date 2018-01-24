function [ solution, residuals, jacobian ] = modelTrackKT_fit( x0, param )

% fit options
options = optimset('Jacobian','on','Display','off','TolX', 1e-10, 'Tolfun', 1e-10,'MaxFunEvals', 1e6, 'MaxIter', 1e6);

shift = [x0(1:7:end-1) x0(2:7:end-1) x0(3:7:end-1)];
scale = x0(4:7:end-1);
angles = [x0(5:7:end-1) x0(6:7:end-1) x0(7:7:end-1)];
% make image parameter matrix
imageParam = zeros(2*param.numPairs,3);
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
    imageParam(idx:idx+1,:) = coords';
end
% make sparse pattern
[x,y,z] = ind2sub(size(param.image),1:numel(param.image));
indx = [x' y' z'];
z = z .* param.psfSigma(1) ./ param.psfSigma(2);
coord = imageParam;
coord(:,3) = coord(:,3) .* param.psfSigma(1) ./ param.psfSigma(2);
sparseP = createSparseDistanceMatrix([x' y' z'],coord,4*param.psfSigma(1));
sparseP = logical(sparseP);
sparseP = sparseP(:,repeatEntries((1:param.numPairs)',7));
sparseP = any(sparseP,2);

% crop image
image = param.image(:);
image = image(sparseP);
indx = indx(sparseP,:);
param.image = image;
param.indx = indx;

% fit
solution = lsqnonlin(@modelTrackKT_fit_objectiveFn,x0,[],[],options,param);

%% obj fn
    function [F,J] = modelTrackKT_fit_objectiveFn(x0, param)
        %MODELTRACKKT_OBJECTIVEFN objective function to be optimised by MODELTRACKKT
        %   EHarry Dec 2012
        % x0 -> [shift; scale; angles; ...; bleachingCorrection]
        bleachingCorrection = x0(end);
        shift = [x0(1:7:end-1) x0(2:7:end-1) x0(3:7:end-1)];
        scale = x0(4:7:end-1);
        angles = [x0(5:7:end-1) x0(6:7:end-1) x0(7:7:end-1)];
        clear x0
        
        % make image parameter matrix
        imageParam = zeros(2*param.numPairs,4);
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
        %         imageParam(:,5) = param.psfSigma(1);
        %         imageParam(:,6) = param.psfSigma(1);
        %         imageParam(:,7) = param.psfSigma(2);
        %         imageParam(:,8) = param.bg;
        
        % fit parameters
        % fitParameters = [{'X1'} {'X2'} {'X3'} {'B'}];
        
        % fit image
        % [~,~,~,~,~,~,resAndGauss,jacobianXYZ] = GaussFitND_noFit(param.image,param.indx,fitParameters,imageParam);
        imageParam = imageParam';
        imageParam = imageParam(:);
        imageParam(end+1) = param.bg;
        [F,jacobianXYZ,gaussIm] = fitNGaussians3D_gaussReturn(imageParam,param.image,param.indx,param.psfSigma);
        jacobianXYZ(:,4:4:end-1) = [];
        
        % set F
        % F = resAndGauss(:,1);
        residuals = F;
        
        % make J
        J = zeros(length(F),(7*param.numPairs)+1);
        for iPair = 1:param.numPairs
            idx = (iPair-1)*7 + 1;
            idxXYZ = (iPair-1)*6 + 1;
            % jacobian rel to shift = sum jabobain x y and z
            J(:,idx:idx+2) = repmat(sum(jacobianXYZ(:,idxXYZ:idxXYZ+5),2),1,3);
            % jacobian rel to scale is (sum jabobain x y and z) .* rotmat * shifted coords
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
            rotCoord = rotationMat * (baseCoords - baseCentre(:,[1 1]));
            J(:,idx+3) = jacobianXYZ(:,idxXYZ).*rotCoord(1,1) + jacobianXYZ(:,idxXYZ+1).*rotCoord(2,1) + jacobianXYZ(:,idxXYZ+2).*rotCoord(3,1) + jacobianXYZ(:,idxXYZ+3).*rotCoord(1,2) + jacobianXYZ(:,idxXYZ+4).*rotCoord(2,2) + jacobianXYZ(:,idxXYZ+5).*rotCoord(3,2);
            % jacobian rel to angles, make dif rot mat and set for each angle
            % wrt phi
            rotationMatPhiD = [0 0 0; 0 -sin(phi) -cos(phi); 0 cos(phi) -sin(phi)];
            rotationMat = rotationMatPsi * rotationMatTheta * rotationMatPhiD;
            rotCoord = rotationMat * (scale(iPair)*eye(3)) * (baseCoords - baseCentre(:,[1 1]));
            J(:,idx+4) = jacobianXYZ(:,idxXYZ).*rotCoord(1,1) + jacobianXYZ(:,idxXYZ+1).*rotCoord(2,1) + jacobianXYZ(:,idxXYZ+2).*rotCoord(3,1) + jacobianXYZ(:,idxXYZ+3).*rotCoord(1,2) + jacobianXYZ(:,idxXYZ+4).*rotCoord(2,2) + jacobianXYZ(:,idxXYZ+5).*rotCoord(3,2);
            % wrt theta
            rotationMatThetaD = [-sin(theta) 0 cos(theta); 0 0 0; -cos(theta) 0 -sin(theta)];
            rotationMat = rotationMatPsi * rotationMatThetaD * rotationMatPhi;
            rotCoord = rotationMat * (scale(iPair)*eye(3)) * (baseCoords - baseCentre(:,[1 1]));
            J(:,idx+5) = jacobianXYZ(:,idxXYZ).*rotCoord(1,1) + jacobianXYZ(:,idxXYZ+1).*rotCoord(2,1) + jacobianXYZ(:,idxXYZ+2).*rotCoord(3,1) + jacobianXYZ(:,idxXYZ+3).*rotCoord(1,2) + jacobianXYZ(:,idxXYZ+4).*rotCoord(2,2) + jacobianXYZ(:,idxXYZ+5).*rotCoord(3,2);
            % wrt psi
            rotationMatPsiD = [-sin(psi) -cos(psi) 0; cos(psi) -sin(psi) 0; 0 0 0];
            rotationMat = rotationMatPsiD * rotationMatTheta * rotationMatPhi;
            rotCoord = rotationMat * (scale(iPair)*eye(3)) * (baseCoords - baseCentre(:,[1 1]));
            J(:,idx+6) = jacobianXYZ(:,idxXYZ).*rotCoord(1,1) + jacobianXYZ(:,idxXYZ+1).*rotCoord(2,1) + jacobianXYZ(:,idxXYZ+2).*rotCoord(3,1) + jacobianXYZ(:,idxXYZ+3).*rotCoord(1,2) + jacobianXYZ(:,idxXYZ+4).*rotCoord(2,2) + jacobianXYZ(:,idxXYZ+5).*rotCoord(3,2);
        end
        %J(:,end) = resAndGauss(:,2) ./ bleachingCorrection;
        J(:,end) = gaussIm ./ bleachingCorrection;
        
        % copy to jacobian
        jacobian = J;
        
        % make sparse pattern
        indx = param.indx;
        indx(:,3) = indx(:,3) .* param.psfSigma(1) ./ param.psfSigma(2);
        sparseP = createSparseDistanceMatrix(indx,coord,3*param.psfSigma(1));
        sparseP = logical(sparseP);
        sparseP = sparseP(:,repeatEntries((1:param.numPairs)',7));
        sparseP(:,end+1) = true;
        
        J = J .* sparseP;
        
    end




end

