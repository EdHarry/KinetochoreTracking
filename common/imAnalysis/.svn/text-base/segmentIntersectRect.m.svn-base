function [xi,yi] = segmentIntersectRect(xS,yS,rect)
%SEGMENTIINTERSECTRECT finds the intersection of a line with a rectangle
%
% SYNOPSIS [xi,yi] = segmentIntersectRect(xS,yS,rect)
%
% INPUT xS : vector with the x-coordinates of the segment
%       yS : vector with the y-coordinates of the segment
%       rect : is a 4-element vector with the form 
%              [xmin ymin width height]. 
%
% OUTPUT xi, yi : x/y coordinate pair of the intersection points
%                 [] is the segment is outside the rectangle
% NOTE  Nonempty vectors are returned even if the segment does
%       not physically intersect the rectangle outline.
%       The term segment is used as 2-point representation of a 
%       line from -infinity to +infinity


% check input data 
if((sum([length(xS),length(yS)]==2)~=2) | length(rect)~=4)
   error('invalid vectors entered');
end;

% get the rectangle corners
C = [rect(1:2);rect(1:2)+[rect(3),0];rect(1:2)+[rect(3),rect(4)];...
      rect(1:2)+[0,rect(4)];rect(1:2)];

xi=[];
yi=[];
for i=1:4
   if(length(xi) < 2)
      % not all the intersection points have been found
      [xAux,yAux] = intersectSegments(C(i:i+1,1),C(i:i+1,2),xS,yS);
      if(isbetween(xAux,yAux,C(i:i+1,1),C(i:i+1,2)))
         xi = [xi xAux];
         yi = [yi yAux];
      end;
   end;
end;


%--------------------------------------------------------------
% local function isbetween

function ret = isbetween(xP,yP,xS,yS)

xmin = min(xS);
xmax = max(xS);
ymin = min(yS);
ymax = max(yS);

ret = (xP >= xmin) & (xP <= xmax) & (yP >= ymin) & (yP <= ymax);