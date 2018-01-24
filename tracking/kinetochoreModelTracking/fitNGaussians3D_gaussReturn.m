function [F,J,gaussIm] = fitNGaussians3D_gaussReturn(x0,image,index,psfSigma)
% Edit of fitNGaussians2D to work in 3D
%   EHarry March 2012

%% ORIGINAL HEADER
% % %FITNGAUSSIANS2D yields F, the difference between an image and a theoretical image produced by N Gaussians, and J, the Jacobian of F.
% % %
% % %SYNOPSIS [F,J] = fitNGaussians2D(x0,image,index,psfSigma)
% % %
% % %INPUT  x0      : initial guess of PSF positions and amplitudes and
% % %                 background noise.
% % %       image   : Image part being analyzed.
% % %       index   : x,y-indices of pixels considered.
% % %       psfSigma: Standard deviation of point spread function (in pixels).
% %
% % %OUTPUT F       : Residuals from fitting an image with supplied
% % %                 Gaussians.
% % %       J       : The Jacobian matrix of F.
% % %       errFlag : 0 if function executes normally, 1 otherwise.
% % %
% % %REMARKS F = model image - real image, important to know if the sign of the
% % %residuals matters.
% % %
% % %Khuloud Jaqaman, August 2005

%% Output

F = [];
J = [];

%% Input

% %check whether correct number of input arguments was used
% if nargin ~= 4
%     disp('--fitNGaussians2D: Incorrect number of input arguments!');
%     return
% end
%check whether correct number of input arguments was used
if nargin ~= 4
    disp('--fitNGaussians3D: Incorrect number of input arguments!');
    return
end

%% Calculating F & J

%extract background intensity from x0 and remove from vector
bgAmp = x0(end);
x0 = x0(1:end-1);

%get number of PSFs considered
% numPSF = length(x0)/3;
numPSF = length(x0)/4;

% %reshape 3nx1 vector x0 into nx3 matrix
% x0 = reshape(x0,3,numPSF);
% x0 = x0';

%reshape 4nx1 vector x0 into nx4 matrix
x0 = reshape(x0,4,numPSF);
x0 = x0';

%extract PSF center positions and amplitudes
% psfPos = x0(:,1:2);
% psfAmp = x0(:,3);

psfPos = x0(:,1:3);
psfAmp = x0(:,4);

%find minimum and maximum pixel indices
minIndxX = min(index(:,1));
maxIndxX = max(index(:,1));
minIndxY = min(index(:,2));
maxIndxY = max(index(:,2));
minIndxZ = min(index(:,3));
maxIndxZ = max(index(:,3));

%determine the contribution of each PSF (assuming amplitude 1) to a
%pixel based on its x-coordinate (needed to calculate F & J)
psfIntegX = zeros(maxIndxX-minIndxX+1,numPSF);
for i=1:numPSF
    psfIntegX(:,i) = GaussListND((minIndxX:maxIndxX)',...
        psfSigma(1),psfPos(i,1));
end

%determine the contribution of each PSF (assuming amplitude 1) to a
%pixel based on its y-coordinate (needed to calculate F & J)
psfIntegY = zeros(maxIndxY-minIndxY+1,numPSF);
for i=1:numPSF
    psfIntegY(:,i) = GaussListND((minIndxY:maxIndxY)',...
        psfSigma(1),psfPos(i,2));
end

%determine the contribution of each PSF (assuming amplitude 1) to a
%pixel based on its z-coordinate (needed to calculate F & J)
psfIntegZ = zeros(maxIndxZ-minIndxZ+1,numPSF);
for i=1:numPSF
    psfIntegZ(:,i) = GaussListND((minIndxZ:maxIndxZ)',...
        psfSigma(2),psfPos(i,3));
end

%get number of pixels in image
numPixel = length(image);

% %get xy-indices relative to minimum
% relIndxX = index(:,1) - minIndxX + 1;
% relIndxY = index(:,2) - minIndxY + 1;
%get xy-indices relative to minimum
relIndxX = index(:,1) - minIndxX + 1;
relIndxY = index(:,2) - minIndxY + 1;
relIndxZ = index(:,3) - minIndxZ + 1;


%calculate the value of F at all pixels
F = (sum(repmat(psfAmp,1,numPixel).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)'.*psfIntegZ(relIndxZ,:)',1))' ...
    + repmat(bgAmp,numPixel,1) - image;
gaussIm = F + image;

%remove pixels with NaN (which means they are out of the cropped image
%area)
indxPixel = find(~isnan(image));

F = F(indxPixel);

if nargout > 1

%calculate the value of each PSF (assuming amplitude 1) at the
%x-coordinates of the corners of all pixels (needed to calculate J)
psfValueX = zeros(maxIndxX-minIndxX+2,numPSF);
for i=1:numPSF
    psfValueX(:,i) = exp(-((minIndxX-0.5:maxIndxX+0.5)'...
        -psfPos(i,1)).^2/2/psfSigma(1)^2);
end

%calculate the value of each PSF (assuming amplitude 1) at the
%y-coordinates of the corners of all pixels (needed to calculate J)
psfValueY = zeros(maxIndxY-minIndxY+2,numPSF);
for i=1:numPSF
    psfValueY(:,i) = exp(-((minIndxY-0.5:maxIndxY+0.5)'...
        -psfPos(i,2)).^2/2/psfSigma(1)^2);
end

%calculate the value of each PSF (assuming amplitude 1) at the
%z-coordinates of the corners of all pixels (needed to calculate J)
psfValueZ = zeros(maxIndxZ-minIndxZ+2,numPSF);
for i=1:numPSF
    psfValueZ(:,i) = exp(-((minIndxZ-0.5:maxIndxZ+0.5)'...
        -psfPos(i,3)).^2/2/psfSigma(2)^2);
end


    %calculate the derivative at all pixels
    %     J = ones(numPixel,3*numPSF+1); %(last column for background amplitude)
    %     J(:,1:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueX(relIndxX,:)-...
    %         psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:); %w.r.t. x
    %     J(:,2:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueY(relIndxY,:)-...
    %         psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:); %w.r.t. y
    %     J(:,3:3:3*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); %w.r.t. amp
    
    J = ones(numPixel,4*numPSF+1); %(last column for background amplitude)
    J(:,1:4:4*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueX(relIndxX,:)-...
        psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:).*psfIntegZ(relIndxZ,:); %w.r.t. x
    J(:,2:4:4*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueY(relIndxY,:)-...
        psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:).*psfIntegZ(relIndxZ,:); %w.r.t. y
    J(:,3:4:4*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueZ(relIndxZ,:)-...
        psfValueZ(relIndxZ+1,:)).*psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); %w.r.t. z
    J(:,4:4:4*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:).*psfIntegZ(relIndxZ,:); %w.r.t. amp
    
    %remove pixels with NaN (which means they are out of the cropped image
    %area)
    J = J(indxPixel,:);
end

%% ~~ the end ~~

%% OLD CODE

% % J = ones(numPixel,3*numPSF+1);
% % F = ones(numPixel,1);
% %
% % for i=1:numPixel %for each pixel
% %
% %     %get xy-indices relative to minimum
% %     relIndxX = index(i,1) - minIndxX + 1;
% %     relIndxY = index(i,2) - minIndxY + 1;
% %
% %     %calculate the value of F
% %     F(i) = sum(psfAmp.*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)') ...
% %         + bgAmp - image(i);
% %
% %     %calculate the derivative wrt x-coordinate
% %     J(i,1:3:3*numPSF) = psfAmp'.*(psfValueX(relIndxX,:)-...
% %         psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:)/psfSigma^2;
% %
% %     %calculate the derivative wrt y-coordinate
% %     J(i,2:3:3*numPSF) = psfAmp'.*(psfValueY(relIndxY,:)-...
% %         psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:)/psfSigma^2;
% %
% %     %calculate the derivative wrt amplitude
% %     J(i,3:3:3*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:);
% %
% %     %since the derivative wrt background intensity = 1, this is already
% %     %accounted for in the initial assignment of J.
% %
% % end



