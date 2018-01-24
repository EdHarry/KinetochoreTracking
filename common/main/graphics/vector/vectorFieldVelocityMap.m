function d=vectorFieldVelocityMap(Mi)
% vectorFieldVelocityMap displays and returns a vector field velocity map
%
% SYNOPSIS   d=vectorFieldVelocityMap(Mi)
%
% INPUT      Mi : vector field M stored in a (nx4)-matrix of the 
%                 form [y0 x0 y x]n, (where (y0,x0) is the base and (y,x) 
%                 is the tip of the vector).
%
% OUTPUT     d  : velocity map
%
% Aaron Ponti 12/06/2002

if nargin~=1
    error('One input parameter expected');
end
if size(Mi,2)~=4
    error('The input parameter Mi must be an (nx4) matrix');
end

% Grid dimensions
ly=length(unique(Mi(:,1)));
lx=length(unique(Mi(:,2)));

% Vectors
Vi=[Mi(:,3)-Mi(:,1) Mi(:,4)-Mi(:,2)];

% Put velocities onto the grid
u=(reshape(Vi(:,1),lx,ly))';
v=(reshape(Vi(:,2),lx,ly))';

% Calculate velocity map
d=sqrt(u.^2+v.^2);

% Display
imshow(d,[]);
colormap('jet');
colorbar
    