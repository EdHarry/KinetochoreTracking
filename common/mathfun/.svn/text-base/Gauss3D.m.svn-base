function out=Gauss3D(x,sigma,fSze);
% Gauss3D	apply a 3 dimensional gauss filter
%
%    out=Gauss3D(x,sigma,fSze);
%
%    INPUT: x      3-dimensional data
%           sigma  of gauss filter [sigmaX sigmaY sigmaZ]
%           fSze   size of the gauss mask [sizeX sizeY sizeZ]
%                  (odd size required for symmetric mask!)
%
%    OUTPUT: out   filtered data

% c: 11/01/00 dT

% fastGauss3D brings equivalent result and is faster
out=fastGauss3D(x,sigma,fSze); % no border correction!

% c=ceil(fSze/2);
% for page=1:fSze(3)
%    for col = 1:fSze(2)
%       for row = 1:fSze(1)
%          mask(row,col,page) = 1/(sqrt(2*pi*dot(sigma,sigma)))^3*...
%             (exp(-1/2*(((c(1)-row)/sigma(1))^2+((c(2)-col)/sigma(2))^2+ ((c(3)-page)/sigma(3))^2)));
%       end;
%    end;
% end;
% % Convolute matrices
% out=convn(x,mask,'same');
% 
% 

