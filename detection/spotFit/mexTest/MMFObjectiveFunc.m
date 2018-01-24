
function [resi,jacobian] = MMFObjectiveFunc(X,maskAmp,maskCoord,sigma)
%MMFOBJECTIVEFUNC is the callback for nonlinear fitter in function mixture
%                 ModelFitting.
%
% SYNOPSIS: [resi,jacobian] = MMFObjectiveFunc(centroid,maskAmp,maskCoord,psfSize)
%
% INPUT centroid : m x n matrix : with m centroid in n-Dimensions.
%       maskAmp  : p x 1 vector with p amplitude corresponding to the
%                  intensity under maskCoord pixel -p.
%       maskCoord : p x n matrix : with m pixel under the mask coordinate in n-Dimensions
%       sigma : n x 1 vector with sigmas used for fitting in n-Dimensions.
%               if sigma is scalar, the same sigma will be used for all
%               dimensions.
% OUTPUT resi : Residual of the mixture gaussian fitting.
%        jacobian : Jacobian of the corresponding fitting matrix.
%
% REMARKS
%
% SEE ALSO mixtureModelFitting, spotMMFit
%
% EXAMPLE
%

% created with MATLAB ver.: 7.14.0.739 (R2012a) on Mac OS X  Version: 10.6.8 Build: 10K549
%
% created by: Jacques Boisvert
% DATE: 28-May-2012
%
% Last revision $Rev: 2723 $ $Date: 2012-05-03 03:33:21 -0400 (Thu, 03 May 2012) $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% reshape X back to 2D
X = reshape(X,size(maskCoord,2),length(X)/size(maskCoord,2))';


%If #Spot > 1, calculate sparse jacobian matrix template
if size(X,1) > 1
    %Input format x1,y1,z1,
    % EHarry, have to warp z coordinate as well
    if length(sigma) == 3
        tmp1 = maskCoord;
        tmp1(:,3) = tmp1(:,3) .* sigma(1) ./ sigma(3);
        tmp2 = X;
        tmp2(:,3) = tmp2(:,3) .* sigma(1) ./ sigma(3);
        sparseJacobian = createSparseDistanceMatrix(tmp1,tmp2,3*sigma(1)) ;
        clear tmp1 tmp2
    else
        sparseJacobian = createSparseDistanceMatrix(maskCoord,X,3*sigma) ;
    end
    %sparseJacobian = sparseJacobian >= 1;
    sparseJacobian = logical(sparseJacobian);
    sparseJacobian = sparseJacobian(:,repeatEntries((1:size(X,1))',size(X,2)));
    %sparseJacobian = full(sparseJacobian);
    %tmpMatrix = zeros(size(sparseJacobian));
    %tmpMatrix(sparseJacobian == 1) = 1;
    %sparseJacobian = tmpMatrix;
else
    sparseJacobian = NaN;
end


%Make sure maskAmp is double.
maskAmp = double(maskAmp);
%% INITIAL ALLOCATION
%To get bg from fitting, fit [G|1] * (a,bg)'
%Ones -> To obtain the bg column from the start. gaussList is 1-n
%gauss + one col of ones for bg.
gaussList = ones(size(maskAmp,1),size(X,1)+1);
%Number of gaussian is equal to the number of spots;
%nGauss = nSpot;
nGauss = size(X,1);
nDim = size(X,2);

%% GAUSS LIST CALCULATION
for gaussIdx = 1:nGauss
    gaussList(:,gaussIdx) = GaussListND(maskCoord,sigma,X(gaussIdx,:));
    %gaussList(:,gaussIdx) = GaussListND_mexCode_mex(maskCoord,sigma,X(gaussIdx,:));
end


%% GAUSS FITTING
%I = (G1,G2,G3,....,1) *(a1;a2;a3;...;bg);
%To fit -> a1;a2..;bg = (G1,G2,G3,...,1) \ I
newAmpBg = gaussList \ maskAmp;

%% RESIDUAL CALCULATION
%Resi : I - I' -> I' = (G1,G2,G3,G4...1) * (a1;a2;a3;...;bg);
resi = maskAmp - ( gaussList * newAmpBg );
bg = newAmpBg(end);
amp = newAmpBg(1:end-1);

%% ASSEMBLY of JACOBIAN MATRIX
%Compute gradiant of function --> aG(x,y,z)+bg = I;
%Preallocation : jacobian is
%[x1,y1,z1,x2,y2,z2...,xn,yn,zn,a1,a2..,an,bg];
%jacobian = zeros(size(gaussList,1), ( nGauss*nDim + nGauss + 1) );
%Increased gaussList
jacobian = repeatEntries(gaussList(:,1:end-1)',nDim)';
%dI/da = G and dI/dbg = 1
jacobian = cat(2,jacobian,gaussList);
%For each coordinate, calculate the first partial derivative--> for
%x,y,z is a*g*(x-c)^2/2*sigma^2 * a*-(x-c)/sigma^2
for cJacobian = 1:(nGauss)
    col = (cJacobian -1) * nDim + 1;%x..y..z
    upperCol = col+nDim-1;
    if(amp(cJacobian) < 0)
        jacobian(:,col:upperCol) = bsxfun(@times,jacobian(:,col:upperCol),-amp(cJacobian));
    else
        jacobian(:,col:upperCol) = bsxfun(@times,jacobian(:,col:upperCol),amp(cJacobian));
    end
    jacobian(:,col:upperCol) = -jacobian(:,col:upperCol) .* bsxfun(@rdivide,bsxfun(@minus,maskCoord,X(cJacobian,:)),(sigma.^2));
end
jacobianMatrix = jacobian;
%Optimisation function only want the jacobian with the x,y,z.
jacobian = jacobian(:,1:nGauss*nDim);

if ~isnan(sparseJacobian)
    jacobian = jacobian .* sparseJacobian;
end

%        fprintf('In objective function\n');
%        fprintf('Centroid, x : %d, y : %d\n',X(1,1),X(1,2));
end

