function [crds,varargout] = plotrect(rect,opt)
%PLOTRECT plots a rectangle
%
% SYNOPSIS:
%    plotrect(rect,opt);
%    crds = plotrect(rect,opt);
%    [crds,objH] = plotrect(rect,opt);
%
% INPUT rect: rectangle [ul(1),ul(2),width,height]
%       opt : plot options (see plot)
%
% OUTPUT:
%    crds: 2x4 matrix with corner coordinates 
%    objH: Handle to the plotted rectangle object.
%
% SEE ALSO plot

if nargout > 2
   error('Too many output arguments. See help plotrect.');
end

x1 = [rect(1),rect(1)+rect(3),rect(1)+rect(3),rect(1),rect(1)];
x2 = [rect(2),rect(2),rect(2)+rect(4),rect(2)+rect(4),rect(2)];

objH = plot(x1,x2,opt);

crds = [x1(1:4);x2(1:4)];

if nargout == 2
   varargout{1} = objH;
end
return;