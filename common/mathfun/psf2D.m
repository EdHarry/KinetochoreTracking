function M = psf2D(pixSize,NA,lambda)

%PSF2D generates a 2D point spread function model (Airy disc)
%
% SYNOPSIS M = psf2D(pixSize,NA,lambda)
%
% INPUT    pixSize : pixel size in object space (in um)
%          NA      : numerical aperture of the objective lens
%          lambda  : wavelength of light (in um)
%
% OUTPUT   M       : filter mask (odd dimensions) representing the Airy disk
%                    the filter is normalized to sum(M(:)) = 1 
%
% NOTE     the function calculates the Airy disk radius according to 
%          R0 = 0.61 * lambda / NA
%          
% REMARK - M is not normalized here!
%
% Alexandre Matov, January 7th, 2003

R0 = 0.61 * lambda / NA;
R1 = 7.02 / 3.83 * R0;    % position of the second root of the Bessel function 

% dim = 2*ceil(R1/pixSize)+1;

[x,y] = meshgrid(-ceil(R1/pixSize):ceil(R1/pixSize));

d = sqrt(x.^2+y.^2)+ 0.01 ;    % shift center by a 1/100 pixel to avoid division-by-0
                               % in the next step
ds = d*pixSize/R0*3.83;
psfs = (besselj(1,ds)./ds);
psf = psfs.*conj(psfs);

% normalization to sum(M(:))=1;
M = psf / sum(psf(:));                                 
%M=psf/max(psf(:));

