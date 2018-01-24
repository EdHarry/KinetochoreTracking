function [Gx,Gy,Gxx,Gyy,LoG,G]=imGaussBench2D(sigma,stretch)
%IMGAUSSBENCH2D returns a bench of 5 spatial filters suggested by Malik for early vision functions
%                 (cf., e.g. Weber & Malik, 1995, IJCV 14:67--81
%
% SYNOPSIS [Gx,Gy,Gxx,Gyy,LoG,G]=imGaussBench2D(sigma,stretch)
%
% INPUT sigma:     filter scale (default 1.0)
%       stretch :  Malik suggests to stretch the directional derivatives 
%                  in the direction perpendicular to the derivative (default 1.4)
%
% OUTPUT Gx, Gy :  directional derivative in x and y
%        Gxx,Gyy:  second directional derivative in x and y
%        LoG    :  Laplacian of Gauss
%        G      :  Gauss kernel itself
%
% SEE ALSO fspecial (Image Processing Toolbox)

% STARTED July-7-1999  G. Danuser

switch nargin
case 0, sigma = 1.0; stretch = 1.4;
case 1, stretch = 1.4;
case 2,
case 3, error('invalid argument list');
end;

if(stretch > 1)
   ssiz=2*ceil(sigma*stretch*5.5)+1;
else
   ssiz=2*ceil(sigma*5.5)+1;
end;
siz = [ssiz ssiz];
[x,y] = meshgrid(-(siz(2)-1)/2:(siz(2)-1)/2,-(siz(1)-1)/2:(siz(1)-1)/2);

% produce Gaussian kernel
C = 1 / (2*pi*sigma^2);
G = C * exp(-(x.*x + y.*y)/(2*sigma*sigma));


C = 1 / (2*pi*sigma^2*stretch);

% produce Gaussian kernels stretched in x and y
hGsx = C*exp(-(x.*x /(sigma*stretch)^2 + y.*y/sigma^2)/2);
hGsy = C*exp(-(x.*x /sigma^2 + y.*y/(sigma*stretch)^2)/2);

% produce directional first derivatives in x and y
Gx = - x.*hGsy / sigma^2;
Gy = - y.*hGsx / sigma^2;

% produce directional second derivatives in x and y
Gxx = -(1-x.^2/sigma^2).*hGsy / sigma^2;
Gyy = -(1-y.^2/sigma^2).*hGsy / sigma^2;

% produce Laplacian of Gauss
LoG = -(1-x.^2/sigma^2).*G / sigma^2 ...
   -(1-y.^2/sigma^2).*G / sigma^2;


