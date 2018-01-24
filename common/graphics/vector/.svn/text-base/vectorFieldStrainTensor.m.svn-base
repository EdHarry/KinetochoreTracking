function [R,S]=vectorFieldStrainTensor(M,Pg,d0,polygon)
% vectorFieldTensor calculates the rotation and strain tensors of an averaged vector field
%
% For a vector field v:  
%    Interpolation with a correlation matrix K : <v>      = K * v      (convolution)
%
%      K = sG*exp(-(dx^2+dy^2)/d0^2) (see vectorFieldInterp)
%
%  dK/dx = -2/D0^2*exp(-(dx^2+dy^2)/d0^2)*dx
%  dK/dy = -2/D0^2*exp(-(dx^2+dy^2)/d0^2)*dy
% 
% where: dx = (xi-xg); dx is a (m x n) matrix
%        dy = (yi-yg); dy is a (m x n) matrix
%        Pi (yi xi)n : coordinated of n vector bases;
%        Pg (yg xg)m : coordinated of m grid points.
%
%            | -Ky*vy  -Kx*vy |
% Tensor T = |                | = 0.5*(T-T') + 0.5*(T+T')
%            | -Ky*vx  -Kx*vx |    rotation     strain
%
% SYNOPSIS   [R,S]=vectorFieldStrainTensor(M,Pg,d0,polygon)
%
% INPUT      M       : vector field, stored in a (nx4)-matrix of the form [y0 x0 y x]n
%                      (where (y0,x0) is the base and (y,x) is the tip of
%                      the vector).
%            Pg      : regular grid points, stored in a (mx2)-matrix of the form [yg xg]m.
%            d0      : parameter for the weight function K=exp(-D.^2/d0),
%                      where D is the distance matrix between all grid
%                      points and all vector (base) positions.
%                      d0 must be a scalar.
%                      Pass d0=0 to let the software set d0 for you (d0 is set equal 3 times
%                      the (larger) grid size).
%                      
%            polygon : (optional - pass polygon=[] to disable). The interpolated vector
%                      can be cropped to remove vectors outside a given region of interest.
%                      To create the polygon use the functions ROIPOLY or
%                      GETLINE. These functions return the polygon vertices
%                      stored in two vectors y and x. Set polygon=[y x] to
%                      use with vectorFieldStrainTensor.
%
% OUTPUT     R       : rotation tensor of the interpolated vector field.
%            S       : strain tensor of the interpolated vector field.
% 
% REMARK: the tensor T (not returned) can be calculated as T=R+S.
%
% DEPENDENCES          vectorFieldStrainTensor uses { createDistanceMatrix (C-MEX function)
%                                                     createDiffMatrix (C-MEX function)
%                                                   } 
%                      vectorFieldStrainTensor is used by { }
%
% Aaron Ponti, 11/18/2002

if nargin~=4
    error('4 input parameters expected.');
end

if prod(size(d0))~=1
    error('The input parameter d0 must be a scalar.');
end

if size(M,2)~=4
    error('The vector field M must be (nx4)');
end

if size(Pg,2)~=2
    error('The grid Pg must be (nx2)');
end

if ~isempty(polygon) 
    if size(polygon,2)~=2
        error('The polygon must be (nx2)');
    end
end

% Vector base positions
Pi=M(:,1:2);

% Vectors
V=[M(:,3)-M(:,1) M(:,4)-M(:,2)];

% Calculate distances
D=createDistanceMatrix(Pg,Pi);

% Calculate d0 if needed
if d0==0
    PgY=unique(Pg(:,1));
    PgX=unique(Pg(:,1));
    d0=3*max(PgY(2)-PgY(1),PgX(2)-PgX(1)); % Set d0 equal 3 times the (larger) grid size
end

% Correlation matrix for the interpolator
G=exp(-D.^2./d0.^2); clear D;

% Calculate weights
sG=sum(G,2); 

% Calculate dy, dx (needed to calculate the derivatives of K)
[dy,dx]=createDiffMatrix(Pg,Pi);

% Calculate dK/dx and dK/dy
Kx=-2./d0.^2.*dx.*G;
Ky=-2./d0.^2.*dy.*G;
clear G dx dy

% Build tensor T (since MATLAB takes the first coordinate as y,
%                         | vy  vx |
% express the tensor as T=|        |.
%                         | uy  ux |
T=zeros(2,2,length(sG));
T(1,1,1:length(sG))=(-Ky*V(:,1))./sG; % vy
T(1,2,1:length(sG))=(-Kx*V(:,1))./sG; % vx
T(2,1,1:length(sG))=(-Ky*V(:,2))./sG; % uy
T(2,2,1:length(sG))=(-Kx*V(:,2))./sG; % ux

% Set the tensor of all vectors outside the passed polygon to [0 0;0 0]
if ~isempty(polygon)
    index=inpolygon(Pg(:,1),Pg(:,2),polygon(:,2),polygon(:,1));
    T(:,:,find(~index),3)=0*T(:,:,find(~index),3);
end

% Rotation tensor
R=0.5*(T-permute(T,[2 1 3])); % 0.5*(T-T')

% Strain tensor
S=0.5*(T+permute(T,[2 1 3])); % 0.5*(T+T')

