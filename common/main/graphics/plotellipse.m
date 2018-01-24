function plotellipse(ctr,r1,r2,rot,opt)
%plotellipse plots an ellipse. 
%
% SYNOPSIS 
%    plotellipse(ctr,r1,r2)
%    plotellipse(ctr,r1,r2,rot)
%    plotellipse(ctr,r1,r2,0,opt)
%    plotellipse(ctr,r1,r2,rot,opt)
%
% INPUT ctr    : Center coordinates [ctr(1),ctr(2)]
%       r1, r2 : Radii of the two axes.
%       rot    : Counter clockwise rotion of the first axis in radian angle.
%                If not rotation, pass 0;
%       opt    : plot options (see plot)
%
% SEE ALSO plot

if nargin < 3
   error('Not enough input arguments.');
end

if nargin > 5
   error('Too many input arguments.');
end

if nargin == 3
   rot = 0;
   opt = [];
end

if nargin == 4
   opt = [];
end

%Sampling angle.
dt = 2*pi/400;
t  = [0:dt:2*pi] + rot;

%X,Y coordinates of sampling points on the ellipse before rotation.
sampleP = [r1*cos(t); r2*sin(t)];

%After rotation.
sampleP = [cos(rot) -sin(rot); sin(rot) cos(rot)]*sampleP;


if isempty(opt)
   plot(sampleP(1,:)+ctr(1),sampleP(2,:)+ctr(2));
else
   plot(sampleP(1,:)+ctr(1),sampleP(2,:)+ctr(2),opt);
end
axis('equal');
