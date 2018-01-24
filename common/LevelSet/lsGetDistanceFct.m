function [dist_fct, dist_fct_matrix] = lsGetDistanceFct(mask_img, grid_coordinates, known_zero_level_points, domain, signed)
% LSGETDISTANCEFCT calculates the distance fct for all gridpoints given zero level points 
%    
%
% SYNOPSIS    [dist_val, dist_val_matrix] = lsGetDistanceFct(mask_img, grid_coordinates, known_zero_level_points, domain, signed)
%
% INPUT         
%               
%               
%               
% 
% OUTPUT        
%
%               
%                           
% DEPENDENCES   lsIntersectSplineLine uses {   fminbnd                             
%                                       }
%
%               lsIntersectSplineLine is used by { lsGetDistanceFct
%                                           }
%
% Matthias Machacek 06/09/04

% calculate the minimal distance for each grid point
h_waitbar = waitbar(0,'Processing');
num_grid_coordinates = size(grid_coordinates,1);
num_p = size(known_zero_level_points,1);
dist_fct = zeros(num_grid_coordinates,1);
for j = 1:num_grid_coordinates
    if mod(j,300) == 0
        waitbar(j/num_grid_coordinates, h_waitbar, num2str(j));
    end
    
    dist = zeros(num_p,1);
    for i = 1:num_p
        dist(i) = sqrt((grid_coordinates(j,1) - known_zero_level_points(i,1))^2 +...
                       (grid_coordinates(j,2) - known_zero_level_points(i,2))^2);
    end
    if mask_img(grid_coordinates(j,2), grid_coordinates(j,1)) > 0 && signed
       level_set_sign = -1; 
    else
       level_set_sign = 1;
    end
    dist_fct(j) = level_set_sign * min(dist);
end
close(h_waitbar);

% put the minimal distances from vetor into matrix form
dist_fct_matrix = reshape(dist_fct,size(domain.y_grid_lines,2),size(domain.x_grid_lines,2));




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