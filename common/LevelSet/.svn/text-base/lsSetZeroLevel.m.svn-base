function [val, val_matrix] = lsSetZeroLevel(grid_coordinates, x_grid_lines, y_grid_lines, x_spline, y_spline, x_size, y_size, x_spacing, y_spacing)

% this function set the zero level curve, interpolation onto the 
% grid points based on the 2-d spline representing initial solution


% discretize edge by intersecting it with the grid lines
% x grid lines 
%option=optimset('display','iter','Diagnostics','on');

x1(1:length(x_grid_lines)) = 0;
x2(1:length(x_grid_lines)) = x_spline.knots(end);
x_intersection_parameter = lsIntersectSplineLine(x_grid_lines, x_spline, x1, x2);

y1(1:length(y_grid_lines)) = 0;
y2(1:length(y_grid_lines)) = y_spline.knots(end);
y_intersection_parameter = lsIntersectSplineLine(y_grid_lines, y_spline, y1, y2);

x_x_intersection = fnval(x_spline,x_intersection_parameter); 
y_x_intersection = fnval(y_spline,x_intersection_parameter);

x_y_intersection = fnval(x_spline,y_intersection_parameter); 
y_y_intersection = fnval(y_spline,y_intersection_parameter);

% all intersection points
x_intersection = cat(2, x_x_intersection, x_y_intersection);
y_intersection = cat(2, y_x_intersection, y_y_intersection);

figure,plot(x_x_intersection,y_x_intersection,'.');
hold on
plot(x_y_intersection,y_y_intersection,'.g');
p=1:x_spline.knots(end);
x_spline = fnval(x_spline, p);
y_spline = fnval(y_spline, p);
plot(x_spline, y_spline,'r');

h_waitbar = waitbar(0,'Processing');
num_grid_coordinates = size(grid_coordinates,1);
for j = 1:num_grid_coordinates
    waitbar(j/num_grid_coordinates, h_waitbar, num2str(j));
    for i = 1:length(x_intersection)
        dist(i) = sqrt((grid_coordinates(j,1) - x_intersection(i))^2 + (grid_coordinates(j,2) - y_intersection(i))^2);
    end
    val(j) = min(dist);
end


val_matrix = reshape(val,size(x_grid_lines,2),size(y_grid_lines,2));
%surface(grid_coordinates(:,1),grid_coordinates(:,2),z);

%D   =  createDistanceMatrix(grid_coordinates(1,:), [x_intersection', y_intersection']);
 
% M and N are the matrices containing the set of point coordinates.
% M and N can represent point positions in 1, 2 and 3D, as follows.
%   
% 
% M=[ y1 x1     and   N=[ y1 x1
% y2 x2              y2 x2
% ...                ...
% ym xm ]            yn xn ]
% OUTPUT   D : distance matrix D=(dij), i=1..m, j=1..n



function f = FunSplineIntersection(x, line, spline)
f = (line - fnval(spline,x))^2;