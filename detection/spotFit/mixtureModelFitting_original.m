function [newCoord,resi,amp,bg,jacobianMatrix] = mixtureModelFitting_original(centroid,maskAmp,maskCoord,sigma,nMaxIter)
%MIXTUREMODELFITTING optimize the localization of the centroid with lsqnonlin
%
% SYNOPSIS: [newCoord,resi,amp,bg] = mixtureModelFitting(coord,amp,bg)
%
% INPUT centroid : m x n matrix : with m center in n-Dimensions.
%       maskAmp  : p x 1 vector with p amplitude corresponding to the
%                  intensity under maskCoord pixel -p.
%       maskCoord : p x n matrix : with m pixel under the mask coordinate in n-Dimensions
%       sigma : n x 1 vector with sigmas used for fitting in n-Dimensions.
%               if sigma is scalar, the same sigma will be used for all
%               dimensions.%
%
%      nMaxIter(opt) : Maximum number of iteration for optimization of
%                      fitting.
%                      def : 0
%
% OUTPUT newCoord : m x n matrix with the optimized centroid
%        resi     : last mixture model fit residual
%        amp      : new estimated intensities
%        jacobianMatrix : Full jacobian matrix from last mixture model
%                         fit. Jacobian has this form : [x1,y1,z1..,zn,a1,a1..,an,bg];
%
% REMARKS
%
% SEE ALSO spotMMFIT.m, MMFObjectiveFunc(nestedFunction)
%
% EXAMPLE
%

% created with MATLAB ver.: 7.14.0.739 (R2012a) on Mac OS X  Version: 10.6.8 Build: 10K549
%
% created by: Jacques Boisvert
% DATE: 24-May-2012
%
% Last revision $Rev: 2723 $ $Date: 2012-05-03 03:33:21 -0400 (Thu, 03 May 2012) $
% Edit by EHarry 19/11/12, make the jacobian coming out of the objective
% function sparse. Forget the JacobPattern option
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ARGUMENT VALIDATION
st = dbstack;
caller = st(1).file;
if nargin < 1 || size(nargin,2) > 3
    error('%s : Please input n x m centroid matrix', caller);
end

if nargin < 2 || size(nargin,2) > 1
    error('%s : Please input p x 1 amplitude list',caller);
end

if nargin < 3 || size(nargin,2) > 3
    error('%s : Please input p x m coordinates corresponding to the amplitude list',caller);
end

if nargin < 4
    error('%s : Please input sigma',caller);
end
if nargin < 5
    nMaxIter = 0;
end


%If #Spot > 1, calculate sparse jacobian matrix template
if size(centroid,1) > 1
   %Input format x1,y1,z1,
   % EHarry, have to warp z coordinate as well
   
   sparseJacobian = createSparseDistanceMatrix(maskCoord,centroid,3*sigma) ;
   sparseJacobian = sparseJacobian >= 1;
   sparseJacobian = sparseJacobian(:,repeatEntries(1:size(centroid,1),size(centroid,1)));
   sparseJacobian = full(sparseJacobian);
   tmpMatrix = zeros(size(sparseJacobian));
   tmpMatrix(sparseJacobian == 1) = 1;
   sparseJacobian = tmpMatrix;
else 
    sparseJacobian = NaN;
end
%% Lsqnonlin option
%Set no bound for now;
lb = [];
ub = [];
if nMaxIter
    if size(centroid,1) > 1
        options = optimset('Jacobian','on','Display','off','TolFun',1e-20,'TolX',1e-3,'MaxIter',nMaxIter,'JacobPattern',sparseJacobian);
    else
        options = optimset('Jacobian','on','Display','off','TolFun',1e-20,'TolX',1e-3,'MaxIter',nMaxIter);
    end
else
    if size(centroid,1) > 1
        options = optimset('Jacobian','on','Display','off','TolFun',1e-20,'TolX',1e-3,'JacobPattern',sparseJacobian);
    else
        options = optimset('Jacobian','on','Display','off','TolFun',1e-20,'TolX',1e-3);
    end
end

%% LSQNONLIN on fit
[newCoord,~,resi] = lsqnonlin(@MMFObjectiveFunc,centroid,lb,ub,options,maskAmp,maskCoord,sigma);

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
        
        %Reshape the centroid to a list [x1,y1,z1;x2,y2,z2...]
%         centroidList = [];
%         for row = 1:nSpot
%             centroidList = cat(1,centroid(1,i*nDim-1:i*nDim));
%         end
%         centroid = centroidList;
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
        
        %        fprintf('In objective function\n');
        %        fprintf('Centroid, x : %d, y : %d\n',X(1,1),X(1,2));
    end

end