function [x,y] = gpts(n)
%GPTS collect a number of points based on mouse clicks
%
% SYNOPSIS [x,y] = gpts(n)
%
% INPUT  n : (optional) collects n points
%            if not set the function gathers an unlimited 
%            number of points until the enter button us hit
%
% OUTPUT coordinates of the points
%
% SEE ALSO ginput()

but = 1;
nPts = 0;

hold on;

if(nargin == 0)
   n = inf;
end;

while(~isempty(but) & (nPts < n))
   [auxX,auxY,but] = ginput(1);
   plot(auxX,auxY,'b.');
   if(~isempty(but))
      nPts = nPts+1;
      x(nPts) = auxX;
      y(nPts) = auxY;
   end;  
end;

hold off;
