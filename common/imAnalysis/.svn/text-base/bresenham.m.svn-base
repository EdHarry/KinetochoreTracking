function xL=bresenham(xS,xE)
%BRESENHAM computes the integer positions on a line between the
% positions xS and xE
%
% SYNOPSIS xL=bresenham(xS,xE)
%
% INPUT xS : coordinate of line start point
%       xE : coordinate of line end point
% 
% OUTPUT xL : 2xn matrix with the coordinates of all the 
%             integer positions on the line 
%
% NOTE : calls tbx_getline in the C toolbox library 

xS = round(xS);
xE = round(xE);

xL = mexGetLine(xS,xE);


