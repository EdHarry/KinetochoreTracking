function PDMap=createDiffPolyMap(I,I2,M,divM,dL,imgSize,polygon)
%
%
%
%
%
% dL reflects the ration    mean of the vector lengths AFTER interpolation
%                           -----------------------------------------------
%                           mean of the vector lengths BEFORE interpolation
%

I=gauss2d(I,5);
I2=gauss2d(I2,5);

% To avoid changes in image intensities due to imaging/illumination
% artifacts, bring the mean intensity values of the two images at the same
% level
Im=mean(I(:)); 
I2m=mean(I2(:));
I2=I2*Im/I2m;

% Calculate vectors
V=[M(:,3)-M(:,1) M(:,4)-M(:,2)];

% Bring the mean vector length back to the value before interpolation
V=V./dL;

% Grid points and dimensions
Pi=divM(:,1:2);
gY=unique(Pi(:,1));
gX=unique(Pi(:,2));
y=length(gY);
x=length(gX);

% Bring the intensity of the first image at 0 level
I0=I-min(I(:));

% Calculate gradient of I0 in x and y
[I0x I0y]=gradient(I0);

% Extract intensities at the grid positions
I0=I0(gY,gX);
I0x=I0y(gY,gX);
I0y=I0y(gY,gX);
I=I(gY,gX);
I2=I2(gY,gX);

% Create (and reshape) Vy and Vx
Vy=(reshape(V(:,1),x,y))';
Vx=(reshape(V(:,2),x,y))';

% Reshape divM
divM=(reshape(divM(:,3),x,y))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate terms of the continuity equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Intensity change
A=I2-I;

% Meshwork compression/dilation
B=I0.*divM;

% Translocation
C=Vy.*I0y+Vx.*I0x;

% Calculate overall map
PDMap=A+B+C;

mx=abs(max([max(A(:)) max(B(:)) max(C(:))]));
% TEMP - plot
subplot(2,2,1)
surf(A);
axis ij
% axis([1 x 1 y -mx mx]);
view(13,68)
title('Intensity change [ A=img2-img ]');
subplot(2,2,2)
surf(B);
axis ij
% axis([1 x 1 y -mx mx]);
view(13,68)
title('Meshwork compression/dilation [ B=I0.*divM ]');
subplot(2,2,3)
surf(C);
axis ij
% axis([1 x 1 y -mx mx]);
view(13,68)
title('Translocation [ C=Vy.*I0y+Vx.*I0x ]');
subplot(2,2,4)
surf(PDMap);
axis ij
% axis([1 x 1 y -mx mx]);
view(13,68)
title('PDMap [ A+B+C ]');
