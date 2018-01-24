function stack = hanningStack(fg,bg,r0,ctr,dim,nse,ecc,nOri)
%HANNINGSTACK generates of a stack of Hanning windows rotating 
% around the center
%
% for further information: > help hanning
%
% SYNOPSIS stack = hanningStack(fg,bg,r0,ctr,dim,ecc,nOri)
%
% INPUT    fg : foreground intensity
%          bg : background intensity
%          r0 : radius of hanning window 
%               (in case of elliptic window r0 = a )
%          ctr: coordinates of the center
%               ctr(1) x- (horizontal) center coordinate
%               ctr(2) y- (vertical) center coordinate
%          dim: image dimension [height, width]
%               CAUTION: dim is compatible with the size() function
%                        but is NOT expressed in the left-handed 
%                        system.
%          nse: variance of the added zero mean Gaussian noise
%          ecc: eccentricity, where ecc = 1 - b^2/a^2 ; 0 < ecc < 1;
%          nOri:number of ellipse orientations -> step = pi/nOri
%
% OUTPUT   stack : stack of frames
%
% SEE ALSO hanning()

dOri = pi/nOri;

stack = hanning(fg,bg,r0,ctr,dim,nse,[ecc, (iOri-1)*dOri]);

for iOri=2:nOri
   I = cat(3,stack, ...
      hanning(fg,bg,r0,ctr,dim,nse,[ecc, (iOri-1)*dOri]));
end;
