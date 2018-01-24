function IG=Gauss2D1(I,SIG);

% Gauss2D1	 applies a two dimensional Gaussian filter
%
% SYNOPSIS   IG=Gauss2d1(I,SIG);
%
% INPUT      I          :   image
%            SIG        :   sigma of the GK
%
% OUTPUT     IG         :   filtered image
%
% REMARKS       difference with function Gauss2D - the Gaussian mask
%               is not mormalized (so the sum is not equal to one)
%
%
% DEPENDENCES   Gauss2D1 uses { filter2 }
%               Gauss2D1 is used by { fsmPrepSubstructMaxima }
%
% Alexandre Matov, November 7th, 2002


R=8; % cutoff radius of the GK

for i = -R:R,
   for j = -R:R,
      P=(i*i+j*j)/2/SIG/SIG;
      M(i+R+1,j+R+1) = exp(-P);
   end
end

% Convolve matrices
IG=filter2(M,I);
