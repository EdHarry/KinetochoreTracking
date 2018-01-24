function intersection_parameter = lsIntersectSplineLine(lines, spline, x1, x2)
% LSINTERSECTSPLINELINE finds all intersections between a set of lines and a spline 
%    
%
%
% SYNOPSIS      intersection_parameter = lsIntersectSplineLine(lines, spline, x1, x2)
%
% INPUT         lines  :  array with line coordinate
%               spline :  one component spline
%               x1     :  lower spline parameter boundary
%               x2     :  upper spline parameter boundary
% 
% OUTPUT        intersection_parameter : parameter of the pline at the
%                                        intersecions
%                           
% DEPENDENCES   lsIntersectSplineLine uses {   fminbnd                             
%                                       }
%
%               lsIntersectSplineLine is used by { lsGetDistanceFct
%                                           }
%
% Matthias Machacek 06/09/04

option=[];
for i = 1:length(lines)
    if x1(i) < x2(i)
        [spline_point(i), fval(i), exitflag(i)]  = fminbnd(@FunSplineIntersection, x1(i), x2(i), option, lines(i), spline);
    else
        fval(i) = 1;
    end
end

if exist('fval','var')
    [index] = find(fval < 0.001);
    intersection_parameter = spline_point(index);
else
    intersection_parameter = [];
end


% lines_red = lines(index);
% x1_red = x1(index);
% x2_red = x2(index);
% if length(lines_red) > 0
%     intersection_parameter1 = lsIntersectSplineLine(lines_red, spline, intersection_parameter + 0.01, x2_red);
% end
% 
% if length(lines_red) > 0
%     intersection_parameter2 = lsIntersectSplineLine(lines_red, spline, x1_red, intersection_parameter - 0.01);
% end
% 
% if exist('intersection_parameter1') & exist('intersection_parameter2')
%     intersection_parameter = cat(2,intersection_parameter, intersection_parameter1, intersection_parameter2);
% elseif exist('intersection_parameter1') 
%     intersection_parameter = cat(2,intersection_parameter, intersection_parameter1);
% elseif exist('intersection_parameter2')
%     intersection_parameter = cat(2,intersection_parameter, intersection_parameter2);   
% end

function f = FunSplineIntersection(x, line, spline)
f = (line - fnval(spline,x))^2;