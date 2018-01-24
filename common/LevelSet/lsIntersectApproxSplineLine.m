function [s_approx, lines_red, ds] = lsIntersectApproxSplineLine(lines, spline, x1, x2)
% LSINTERSECTSPLINELINE finds all approximately intersections between set of lines and spline 
%    
%
%
% SYNOPSIS      [s_approx, lines_red] = lsIntersectApproxSplineLine(lines, spline, x1, x2)
%
% INPUT         lines  :  array with line coordinate
%               spline :  one component spline
%               x1     :  lower spline parameter boundary
%               x2     :  upper spline parameter boundary
% 
% OUTPUT        s_approx: parameter of the pline at the intersecions
%
%               lines_red: the lines intersecting the curve
%                           
% DEPENDENCES   lsIntersectSplineLine uses {   fminbnd                             
%                                       }
%
%               lsIntersectSplineLine is used by { lsGetDistanceFct
%                                           }
%
% Matthias Machacek 06/09/04


% find the approximate intersection

s_approx = [];
lines_red = [];
for i = 1:length(lines)
   ds = (x2(i) -  x1(i)) / 1000;
   s = x1(i):ds:x2(i);
   f = (lines(i) - fnval(spline,s)).^2;
   local_min = imregionalmin(f);
   [local_min_dist, local_min_index] = find(local_min);
   [min_dist, min_index] = find(f(local_min_index) <= 1); 
   s_approx = cat(2,s_approx,s(local_min_index(min_index)));
   if length(min_index) > 0
      for j=1:length(min_index)
         lines_red(end+1) = lines(i);
      end
   end
end



