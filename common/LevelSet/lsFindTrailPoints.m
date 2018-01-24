function trialGridPoints = lsFindTrailPoints(x_X_intersection, y_X_intersection, x_Y_intersection, y_Y_intersection, domain)
% LSFINDGRIDPOINTS it finds the trail points on a grid given a boundary
% 
%
% SYNOPSIS      trialGridPoints = lsFindTrailPoints(x_grid, y_grid, intersections) 
%
% INPUT         x_X_intersection    :
%               y_X_intersection    :
%               x_Y_intersection    :
%               y_Y_intersection    :
%               domain              :         
% 
% OUTPUT        trialGridPoints       
%                     
%                           
% DEPENDENCES    lsLineMatching uses { remove_multiple_entries
%                                    } 
%
%                lsLineMatching is used by {
%                                             }
%
%
% Matthias Machacek 6/15/04

% get the trail points on the x-grid lines
num_X_intersections = length(x_X_intersection);
num_y_grid_lines    = length(domain.y_grid_lines);

for i=1:num_X_intersections
    j=1;
    while domain.y_grid_lines(j) < y_X_intersection(i) & j <= num_y_grid_lines 
        j=j+1;
    end
    %if j < num_y_grid_lines
%         trialGridPoints_X(2*i-1,2)    = domain.y_grid_lines(j-1);
%         trialGridPoints_X(2*i-1,1)    = x_X_intersection(i);
%         trialGridPoints_X(2*i,2)      = domain.y_grid_lines(j);
%         trialGridPoints_X(2*i,1)      = x_X_intersection(i);    
        trialGridPoints_X(2*i-1,2)    = j-1;
        trialGridPoints_X(2*i-1,1)    = find(abs(domain.x_grid_lines - x_X_intersection(i)) < 0.001);
        trialGridPoints_X(2*i,2)      = j;
        trialGridPoints_X(2*i,1)      = trialGridPoints_X(2*i-1,1);  
    %end
end



% get the trail points on the y-grid lines
num_Y_intersections = length(y_Y_intersection);
num_x_grid_lines    = length(domain.x_grid_lines);

for i=1:num_Y_intersections
    j=1;
    while domain.x_grid_lines(j) < x_Y_intersection(i) & j <= num_x_grid_lines 
        j=j+1;
    end
    %if j < num_y_grid_lines
%         trialGridPoints_Y(2*i-1,1)    = domain.x_grid_lines(j-1);
%         trialGridPoints_Y(2*i-1,2)    = y_Y_intersection(i);
%         trialGridPoints_Y(2*i,1)      = domain.x_grid_lines(j);
%         trialGridPoints_Y(2*i,2)      = y_Y_intersection(i);  
        trialGridPoints_Y(2*i-1,1)    = j-1;
        trialGridPoints_Y(2*i-1,2)    = find(abs(domain.y_grid_lines - y_Y_intersection(i)) < 0.001);
        trialGridPoints_Y(2*i,1)      = j;
        trialGridPoints_Y(2*i,2)      = trialGridPoints_Y(2*i-1,2);          
    %end
end

trialGridPoints_temp = cat(1, trialGridPoints_X, trialGridPoints_Y);

% remove double entries
[trialGridPoints(:,1),trialGridPoints(:,2)] =...
    remove_multiple_entries(trialGridPoints_temp(:,1),trialGridPoints_temp(:,2));

figure,plot(domain.x_grid_lines(trialGridPoints(:,1)), domain.y_grid_lines(trialGridPoints(:,2)),'x');
hold on 
plot(x_X_intersection,y_X_intersection,'o');
plot(x_Y_intersection,y_Y_intersection,'ro');





