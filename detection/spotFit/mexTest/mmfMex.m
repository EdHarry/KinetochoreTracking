function [ jacobian,resi,amp,bg,nGauss,nDim ] = mmfMex( maskAmp,maskCoord,X,sigma )
%MMFMEX Summary of this function goes here
%   Detailed explanation goes here

%% INITIAL ALLOCATION
%To get bg from fitting, fit [G|1] * (a,bg)'
%Ones -> To obtain the bg column from the start. gaussList is 1-n
%gauss + one col of ones for bg.
%gaussList = ones(size(maskAmp,1),size(X,1)+1);
%Number of gaussian is equal to the number of spots;
%nGauss = nSpot;
nGauss = size(X,1);
nDim = size(X,2);

%gaussList = GaussListND_mexCode_(maskCoord,sigma,X);

%nGauss = size(centerInput,1);
nCoords = size(maskCoord,1);

gaussListOutput = ones(nCoords,nGauss+1);

sigma2 = repmat(sigma,[nCoords,1]);

tmp = 0.5./sigma2;

normN = (2*pi)^(0.5*nDim)*prod(sigma2(1,:));

sq2 = sqrt(2);

for gaussIdx = 1:nGauss
    
    
    
    
    center2 = repmat(X(gaussIdx,:), [nCoords,1]);
    
    
    
    %======================
    % CALC GAUSSLIST
    %======================
    
    
    % instead of calculating Gauss-values for very complicated geometries, we
    % make a coordinate transformation so that we can use sigma=1 in all
    % dimensions
    
    
    
    % 0.5*erfc(-(x+0.5)/sqrt(2))-0.5*erfc(-(x-0.5)/sqrt(2)) gives the integral on the
    % pixel at 1 of a Gaussian with mean 0 and sigma 1
    
    
    %center3 = center2(:,1);
    %clear center2
    
    % convert coordList to 0/1
    %coordList2 = (coordList(1:nCoords,1:nDims) - center2(1:nCoords,1:nDims))./sigma2(1:nCoords,1:nDims);
    coordList2 = (maskCoord- center2)./sigma2;
    %clear coordList center3
    
    % double coordList as preparation for erfc
    %fixed bug: must divide the 0.5 by sigma - KJ
    %coordList2 = cat(3,coordList2-0.5./sigma2(1:nCoords,1:nDims), coordList2+0.5./sigma2(1:nCoords,1:nDims));
    
    coordList2 = cat(3,coordList2-tmp, coordList2+tmp);
    
    % calculate gaussList
    %Jonas was missing the minus sign in erfc. I corrected that - KJ
    gaussList = diff(0.5 * erfc(-coordList2/sq2),1,3);
    gaussList = prod(gaussList,2);
    
    % norm gaussList
    gaussListOutput(:,gaussIdx) = gaussList(:,1) * normN;
    
end


%% GAUSS LIST CALCULATION
% for gaussIdx = 1:nGauss
%     tmp = GaussListND_mexCode(maskCoord,sigma,X(gaussIdx,:));
%     gaussList(:,gaussIdx) = tmp(:,1);
%     %gaussList(:,gaussIdx) = GaussListND_mexCode_mex(maskCoord,sigma,X(gaussIdx,:));
% end


%% GAUSS FITTING
%I = (G1,G2,G3,....,1) *(a1;a2;a3;...;bg);
%To fit -> a1;a2..;bg = (G1,G2,G3,...,1) \ I
newAmpBg = gaussListOutput \ maskAmp;

%% RESIDUAL CALCULATION
%Resi : I - I' -> I' = (G1,G2,G3,G4...1) * (a1;a2;a3;...;bg);
resi = maskAmp - ( gaussListOutput * newAmpBg );
bg = newAmpBg(end);
amp = newAmpBg(1:end-1);


sigmaSq = (sigma.^2);

%% ASSEMBLY of JACOBIAN MATRIX
%Compute gradiant of function --> aG(x,y,z)+bg = I;
%Preallocation : jacobian is
%[x1,y1,z1,x2,y2,z2...,xn,yn,zn,a1,a2..,an,bg];
%jacobian = zeros(size(gaussList,1), ( nGauss*nDim + nGauss + 1) );
%Increased gaussList
jacobian = repeatEntries_mex(gaussListOutput(:,1:end-1)',nDim)';
%dI/da = G and dI/dbg = 1
jacobian = cat(2,jacobian,gaussListOutput);
%For each coordinate, calculate the first partial derivative--> for
%x,y,z is a*g*(x-c)^2/2*sigma^2 * a*-(x-c)/sigma^2
for cJacobian = 1:(nGauss)
    col = (cJacobian -1) * nDim + 1;%x..y..z
    upperCol = col+nDim-1;
    if(amp(cJacobian) < 0)
        %jacobian(:,col:upperCol) = bsxfun(@times,jacobian(:,col:upperCol),-amp(cJacobian));
        jacobian(:,col:upperCol) = jacobian(:,col:upperCol) * -amp(cJacobian);
    else
        %jacobian(:,col:upperCol) = bsxfun(@times,jacobian(:,col:upperCol),amp(cJacobian));
        jacobian(:,col:upperCol) = jacobian(:,col:upperCol) * amp(cJacobian);
    end
    jacobian(:,col:upperCol) = -jacobian(:,col:upperCol) .* bsxfun(@rdivide,bsxfun(@minus,maskCoord,X(cJacobian,:)),sigmaSq);
    %tmp = -jacobian(:,col:upperCol) .* ((maskCoord - X(cJacobian,:)) ./ sigmaSq);
    %jacobian(:,col:upperCol) = tmp(:,1);
end



end

