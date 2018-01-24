function [newCoord,resi,amp,bg,jacobianMatrix] = mixtureModelFitting_kjMethod(centroid,maskAmp,maskCoord,sigma,nMaxIter)
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
% new, eharry will change the way the jacobian is calculated to make it
% idential with finNGaussians.m
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
    if length(sigma) == 3
        tmp1 = maskCoord;
        tmp1(:,3) = tmp1(:,3) .* sigma(1) ./ sigma(3);
        tmp2 = centroid;
        tmp2(:,3) = tmp2(:,3) .* sigma(1) ./ sigma(3);
        sparseJacobian = createSparseDistanceMatrix(tmp1,tmp2,3*sigma(1)) ;
        clear tmp1 tmp2
    else
        sparseJacobian = createSparseDistanceMatrix(maskCoord,centroid,3*sigma) ;
    end
    %sparseJacobian = sparseJacobian >= 1;
    sparseJacobian = logical(sparseJacobian);
    sparseJacobian = sparseJacobian(:,repeatEntries((1:size(centroid,1))',size(centroid,2)));
    %sparseJacobian = full(sparseJacobian);
    %tmpMatrix = zeros(size(sparseJacobian));
    %tmpMatrix(sparseJacobian == 1) = 1;
    %sparseJacobian = tmpMatrix;
else
    sparseJacobian = NaN;
end

%% Lsqnonlin option
%Set no bound for now;
lb = [];
ub = [];
if nMaxIter
    %     if size(centroid,1) > 1
    %         options = optimset('Jacobian','on','Display','off','TolFun',1e-20,'TolX',1e-3,'MaxIter',nMaxIter,'JacobPattern',sparseJacobian);
    %     else
    %         options = optimset('Jacobian','on','Display','off','TolFun',1e-20,'TolX',1e-3,'MaxIter',nMaxIter);
    %     end
    options = optimset('Jacobian','on','Display','off','TolFun',1e-20,'TolX',1e-8,'MaxIter',nMaxIter);
else
    %     if size(centroid,1) > 1
    %         options = optimset('Jacobian','on','Display','off','TolFun',1e-20,'TolX',1e-3,'JacobPattern',sparseJacobian);
    %     else
    %         options = optimset('Jacobian','on','Display','off','TolFun',1e-20,'TolX',1e-3);
    %     end
    options = optimset('Jacobian','on','Display','off','TolFun',1e-20,'TolX',1e-8);
end

%% LSQNONLIN on fit
nDims = size(centroid,2);
centroid = centroid';
centroid = centroid(:);
[newCoord,~,resi] = lsqnonlin(@MMFObjectiveFunc,centroid,lb,ub,options,maskAmp,maskCoord,sigma,sparseJacobian);
newCoord = reshape(newCoord,nDims,length(newCoord)/nDims)';

    function [resi,jacobian] = MMFObjectiveFunc(X,maskAmp,maskCoord,sigma,sparsePattern)
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
        
        X = reshape(X,size(maskCoord,2),length(X)/size(maskCoord,2))';
        
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
        
        %% GAUSS LIST IN EACH DIMENSION
        %find minimum and maximum pixel indices
        minIdx = min(maskCoord);
        maxIdx = max(maskCoord);
        
        %determine the contribution of each PSF (assuming amplitude 1) to a
        %pixel based on its x-coordinate (needed to calculate F & J)
        psfIntegX = zeros(maxIdx(1)-minIdx(1)+1,nGauss);
        for i=1:nGauss
            psfIntegX(:,i) = GaussListND((minIdx(1):maxIdx(1))',...
                sigma(1),X(i,1));
        end
        
        if length(sigma) == 1
            sigma = [sigma sigma];
        end
        
        %determine the contribution of each PSF (assuming amplitude 1) to a
        %pixel based on its y-coordinate (needed to calculate F & J)
        psfIntegY = zeros(maxIdx(2)-minIdx(2)+1,nGauss);
        for i=1:nGauss
            psfIntegY(:,i) = GaussListND((minIdx(2):maxIdx(2))',...
                sigma(2),X(i,2));
        end
        
        if nDim == 3
            
            if length(sigma) == 2
                sigma = [sigma sigma(1)];
            end
            
            %determine the contribution of each PSF (assuming amplitude 1) to a
            %pixel based on its z-coordinate (needed to calculate F & J)
            psfIntegZ = zeros(maxIdx(3)-minIdx(3)+1,nGauss);
            for i=1:nGauss
                psfIntegZ(:,i) = GaussListND((minIdx(3):maxIdx(3))',...
                    sigma(3),X(i,3));
            end
            
            
        end
        
        % get relative indexes
        relIdx = maskCoord - repmat(minIdx,size(maskCoord,1),1) + 1;
        
        gaussList(:,1:nGauss) = psfIntegX(relIdx(:,1),:) .* psfIntegY(relIdx(:,2),:);
        
        if nDim == 3
            gaussList(:,1:nGauss) = gaussList(:,1:nGauss) .* psfIntegZ(relIdx(:,3),:);
        end
        
        %         %% GAUSS LIST CALCULATION
        %         for gaussIdx = 1:nGauss
        %             gaussList(:,gaussIdx) = GaussListND(maskCoord,sigma,X(gaussIdx,:));
        %         end
        
        
        %% GAUSS FITTING
        %I = (G1,G2,G3,....,1) *(a1;a2;a3;...;bg);
        %To fit -> a1;a2..;bg = (G1,G2,G3,...,1) \ I
        newAmpBg = gaussList \ maskAmp;
        
        %% RESIDUAL CALCULATION
        %Resi : I - I' -> I' = (G1,G2,G3,G4...1) * (a1;a2;a3;...;bg);
        %resi = maskAmp - ( gaussList * newAmpBg );
        resi = ( gaussList * newAmpBg ) - maskAmp;
        bg = newAmpBg(end);
        amp = newAmpBg(1:end-1);
        
        %% ASSEMBLY of JACOBIAN MATRIX
        %         %Compute gradiant of function --> aG(x,y,z)+bg = I;
        %         %Preallocation : jacobian is
        %         %[x1,y1,z1,x2,y2,z2...,xn,yn,zn,a1,a2..,an,bg];
        %         %jacobian = zeros(size(gaussList,1), ( nGauss*nDim + nGauss + 1) );
        %         %Increased gaussList
        %         jacobian = repeatEntries(gaussList(:,1:end-1)',nDim)';
        %         %dI/da = G and dI/dbg = 1
        %         jacobian = cat(2,jacobian,gaussList);
        %         %For each coordinate, calculate the first partial derivative--> for
        %         %x,y,z is a*g*(x-c)^2/2*sigma^2 * a*-(x-c)/sigma^2
        %         for cJacobian = 1:(nGauss)
        %             col = (cJacobian -1) * nDim + 1;%x..y..z
        %             upperCol = col+nDim-1;
        %             if(amp(cJacobian) < 0)
        %                 jacobian(:,col:upperCol) = bsxfun(@times,jacobian(:,col:upperCol),-amp(cJacobian));
        %             else
        %                 jacobian(:,col:upperCol) = bsxfun(@times,jacobian(:,col:upperCol),amp(cJacobian));
        %             end
        %             jacobian(:,col:upperCol) = -jacobian(:,col:upperCol) .* bsxfun(@rdivide,bsxfun(@minus,maskCoord,X(cJacobian,:)),(sigma.^2));
        %         end
        %         jacobianMatrix = jacobian;
        %         %Optimisation function only want the jacobian with the x,y,z.
        %         jacobian = jacobian(:,1:nGauss*nDim);
        
        %calculate the value of each PSF (assuming amplitude 1) at the
        %x-coordinates of the corners of all pixels (needed to calculate J)
        psfValueX = zeros(maxIdx(1)-minIdx(1)+2,nGauss);
        for i=1:nGauss
            psfValueX(:,i) = exp(-((minIdx(1)-0.5:maxIdx(1)+0.5)'...
                -X(i,1)).^2/2/sigma(1)^2);
        end
        
        %calculate the value of each PSF (assuming amplitude 1) at the
        %y-coordinates of the corners of all pixels (needed to calculate J)
        psfValueY = zeros(maxIdx(2)-minIdx(2)+2,nGauss);
        for i=1:nGauss
            psfValueY(:,i) = exp(-((minIdx(2)-0.5:maxIdx(2)+0.5)'...
                -X(i,2)).^2/2/sigma(2)^2);
        end
        
        if nDim == 3
            
            %calculate the value of each PSF (assuming amplitude 1) at the
            %z-coordinates of the corners of all pixels (needed to calculate J)
            psfValueZ = zeros(maxIdx(3)-minIdx(3)+2,nGauss);
            for i=1:nGauss
                psfValueZ(:,i) = exp(-((minIdx(3)-0.5:maxIdx(3)+0.5)'...
                    -X(i,3)).^2/2/sigma(3)^2);
            end
            
        end
        
        jacobian = zeros(size(gaussList,1), ( nGauss*nDim) );
        jacobian(:,1:nDim:end) = repmat(amp',size(gaussList,1),1).*(psfValueX(relIdx(:,1),:)-psfValueX(relIdx(:,1)+1,:)).*psfIntegY(relIdx(:,2),:); %w.r.t. x
        if nDim == 3
            jacobian(:,1:nDim:end) = jacobian(:,1:nDim:end).*psfIntegZ(relIdx(:,3),:);
        end
        jacobian(:,2:nDim:end) = repmat(amp',size(gaussList,1),1).*(psfValueY(relIdx(:,2),:)-psfValueY(relIdx(:,2)+1,:)).*psfIntegX(relIdx(:,1),:); %w.r.t. y
        if nDim == 3
            jacobian(:,2:nDim:end) = jacobian(:,2:nDim:end).*psfIntegZ(relIdx(:,3),:);
            jacobian(:,3:nDim:end) = repmat(amp',size(gaussList,1),1).*(psfValueZ(relIdx(:,3),:)-psfValueZ(relIdx(:,3)+1,:)).*psfIntegX(relIdx(:,1),:).*psfIntegY(relIdx(:,2),:); %w.r.t. y
        end
  
        % add jacobian for apmlitude and background
        jacobianMatrix = cat(2,jacobian,gaussList);
        
        % finally make sparse if needed
        if ~isnan(sparsePattern)
            jacobian = jacobian .* sparsePattern;
        end
        
        %        fprintf('In objective function\n');
        %        fprintf('Centroid, x : %d, y : %d\n',X(1,1),X(1,2));
    end

end