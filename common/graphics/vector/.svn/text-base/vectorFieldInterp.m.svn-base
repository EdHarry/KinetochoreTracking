function Mi=vectorFieldInterp(M,Pg,d0,polygon)
% vectorFieldInterp interpolates a vector field on a user-specified grid
%
% For a vector field v:  
%    Interpolation with a correlation matrix K : <v> = K * v      (convolution)
%    where:  K = sG*exp(-(dx^2+dy^2)/d0^2), sG = weigth vector
%
% SYNOPSIS   Mi=vectorFieldInterp(M,Pg,d0,polygon)
%
% INPUT      M       : vector field, stored in a (nx4)-matrix of the form [y0 x0 y x]n
%                      (where (y0,x0) is the base and (y,x) is the tip of
%                      the vector).
%            Pg      : regular grid points, stored in a (mx2)-matrix of the form [yg xg]m.
%            d0      : parameter for the weight function G=exp(-D.^2/(1+d0^2)),
%                      where D is the distance matrix between all grid
%                      points and all vector (base) positions.
%                      d0 can be a scalar or a vector with size (nx1). See
%                      updateD0fromDiv.
%            polygon : (optional - pass polygon=[] to disable). The interpolated vector
%                      can be cropped to remove vectors outside a given region of interest.
%                      To create the polygon use the functions ROIPOLY or
%                      GETLINE. These functions return the polygon vertices
%                      stored in two vectors y and x. Set polygon=[x y] to
%                      use with vectorFieldInterp.
%
% OUTPUT     Mi      : interpolated vector field.
%
% DEPENDENCES          vectorFieldInterp uses { createDistanceMatrix (C-MEX function) }
%                      vectorFieldInterp is used by { }
%
% Aaron Ponti, 11/18/2002

% Check input
if nargin~=4
    error('4 input parameters expected.');
end

if prod(size(d0))~=1 & prod(size(d0))~=size(Pg,1)
    error('The input parameter d0 must be a scalar or a (nx1) vector.');
end

if size(M,2)~=4
    error('The vector field M must be (nx4)');
end

if size(Pg,2)~=2
    error('The grid Pg must be (nx2)');
end

% Vector base positions
Pi=M(:,1:2);

% Vectors
V=[M(:,3)-M(:,1) M(:,4)-M(:,2)];

% Calculate distances
D=createDistanceMatrix(Pg,Pi);

% Correlation matrix (d0 may be a scalar or a vector)
G=zeros(size(D));
if size(d0)==[1 1]
    G=exp(-D.^2./d0^2);
else
    for i=1:size(D,1)  
        G(i,:)=exp(-D(i,:).^2./d0(i,1).^2);
    end
end
clear D % Not needed any more

% Interpolate
Vi=[G*V(:,1) G*V(:,2)];

% Normalize
sG=sum(G,2);
sG(find(sG==0))=1; % Prevent division by zero
Vi=[Vi(:,1)./sG Vi(:,2)./sG];

% Mi is the interpolated M
Mi=[Pg Pg+Vi];

% Set all vectors outside the passed polygon to 0
if ~isempty(polygon)
    index=inpolygon(Mi(:,1),Mi(:,2),polygon(:,2),polygon(:,1));
    Mi(find(~index),3)=Mi(find(~index),1);
    Mi(find(~index),4)=Mi(find(~index),2);
end
