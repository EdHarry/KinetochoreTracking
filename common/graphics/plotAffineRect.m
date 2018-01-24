function [crds] = plotAffineRect(ctr,dim,a,opts)
%PLOTAFFINERECT plots a rectangle which is deformed with an affine transform
%
% SYNOPSIS plotAffineRect(ctr,dim,a,opts)
%
% INPUT ctr : center of the rectangle
%       dim : dimension [width, height] of the rectangle
%       a   : affine matrix
%       opts: plot options (see also plot)
%
% OUTPUT crds : 2x4 matrix with coordinates of the corners
%
% SEE ALSO plot

hDim = dim/2;

% upper left coordinate
auxDim = [-hDim(1);-hDim(2)];
aux = a * auxDim;
ul = ctr + aux';
%upper right
auxDim(1) = hDim(1);
aux = a * auxDim;
ur = ctr + aux';
%lower right
auxDim(2) = hDim(2);
aux = a * auxDim;
lr = ctr + aux';
%lower left 
auxDim(1) = -hDim(1);
aux = a * auxDim;
ll = ctr + aux';

x1 = [ul(1),ur(1),lr(1),ll(1),ul(1)];
x2 = [ul(2),ur(2),lr(2),ll(2),ul(2)];
plot(x1,x2,opts);

crds = [x1(1:4);x2(1:4)];

