function plotmask(mask,ctr,opt)
%PLOTMASK plots a mask
%
% SYNOPSIS plotrect(mask,ctr,opt)
%
% INPUT mask: binary mask 
%       ctr : center of the mask
%       opt : plot options (see plot); must be 'c.', where c is 
%             a MATLAB standard color definition
%
% SEE ALSO plot

sze = size(mask);
[xi,yi] = bwExtractCrds(mask);
lh=plot(xi-1+ctr(1)-(sze(2)-1)/2,yi-1+ctr(2)-(sze(1)-1)/2,opt);
set(lh,'MarkerSize',get(lh,'MarkerSize')/3);
