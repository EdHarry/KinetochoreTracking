function dist_fct = lsGetDistanceFctVec(mask_img, grid_coordinates, known_zero_level_points, domain, signed)
% LSGETDISTANCEFCT calculates the distance fct for all gridpoints given zero level points 
%    
%
% SYNOPSIS   dist_fct  = lsGetDistanceFct(mask_img, grid_coordinates, known_zero_level_points, domain, signed)
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
for j = 1:num_grid_coordinates
    waitbar(j/num_grid_coordinates, h_waitbar, num2str(j));
    for i = 1:size(known_zero_level_points,1)
        dist(i) = sqrt((domain.x_grid_lines(grid_coordinates(j,1)) - known_zero_level_points(i,1))^2 +...
                       (domain.y_grid_lines(grid_coordinates(j,2)) - known_zero_level_points(i,2))^2);
    end
    if mask_img(grid_coordinates(j,2), grid_coordinates(j,1)) > 0 & signed
       level_set_sign = -1; 
    else
       level_set_sign = 1;
    end
    dist_fct(j) = level_set_sign * min(dist);
end
