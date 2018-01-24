function dist_fct_matrix = lsFastMarching(dist_fct_matrix, trial_grid_points, domain, LARGE_NUMBER)
% LSFASTMARCHING solves |grad F|=1 on a grid given "trial_grid_points" 
%
%
%
%      Input:  dist_fct_matrix:     memory for solution F  
%              trial_grid_points:   initial values for of "dist_fct_matrix"
%              domain:              structure to describe computational grid
%              LARGE_NUMBER:        large number (1000000)
%
%      Output: dist_fct_matrix:     solution F (distance transformation)
%
%
%   Matthias Machacek Oct. 25 2004


h = domain.x_spacing;

num_x_grid_lines = length(domain.x_grid_lines);
num_y_grid_lines = length(domain.y_grid_lines);

% number of trail points 
num_trail = size(trial_grid_points, 1);

% status matrix
% 0 unknown point
% 1 trail point 
% 2 known point
status_matrix = zeros(num_x_grid_lines,num_y_grid_lines); 

for i=1:size(trial_grid_points,1)
    status_matrix(trial_grid_points(i,1), trial_grid_points(i,2)) = 1; 
    dist_fct_matrix(trial_grid_points(i,1), trial_grid_points(i,2)) =...
                                            trial_grid_points(i,3); 
end


figure
while num_trail > 0
    % find trail point with smallest value
    [min_val, min_index] = min(trial_grid_points(:,3));
    i_min = trial_grid_points(min_index,1);
    j_min = trial_grid_points(min_index,2);


    % remove this point from the trail point vector
    trial_grid_points(min_index,:) = [];
    % tag this point as known
    status_matrix(i_min, j_min) = 2;
    num_trail = num_trail-1;
    % if the neigbourhing points are not known points tag
    % them as trail points and add them to the trail list
    % a: right
    % b: left
    % c: upper
    % d: lower
    
    % update the value of the right neighbour
    if i_min < num_x_grid_lines
        if status_matrix(i_min+1, j_min) ~=2
            if status_matrix(i_min+1, j_min) == 0
                % this point is not yet in the trail list, so put it there
                num_trail = num_trail+1;
                status_matrix(i_min+1, j_min) = 1;

                trial_grid_points(num_trail,1) = i_min+1;
                trial_grid_points(num_trail,2) = j_min;
            end

            % get the new dist fct value of this trail point
            if i_min < num_x_grid_lines-1
                a = dist_fct_matrix(i_min+2  ,j_min);
            else
                a = LARGE_NUMBER;
            end
            b = dist_fct_matrix(i_min, j_min);
            if j_min < num_y_grid_lines
                c = dist_fct_matrix(i_min+1, j_min+1);
            else 
                c = LARGE_NUMBER;
            end
            if j_min > 1
                d = dist_fct_matrix(i_min+1, j_min-1);
            else
                d = LARGE_NUMBER;
            end
            x0 = dist_fct_matrix(i_min,j_min);
            trial_grid_points(num_trail,3) = fzero(@eikonal,x0,[],a,b,c,d,h);
            dist_fct_matrix(i_min+1, j_min) = trial_grid_points(num_trail,3);
            
            % update the extension velocity by solving  |grad F grad phi| = 0
            
            
        end
    end
    % update the value of the left neighbour
    if i_min > 1
        if status_matrix(i_min-1, j_min) ~=2
            if status_matrix(i_min-1, j_min) == 0
                % this point is not yet in the trail list, so put it there
                num_trail = num_trail+1;
                status_matrix(i_min-1, j_min) = 1;

                trial_grid_points(num_trail,1) = i_min-1;
                trial_grid_points(num_trail,2) = j_min;
            end
            
            % get the new dist fct value of this trail point
            a = dist_fct_matrix(i_min, j_min);
            if i_min > 2
                b = dist_fct_matrix(i_min-2, j_min);
            else
                b = LARGE_NUMBER;
            end
            if j_min < num_y_grid_lines
                c = dist_fct_matrix(i_min-1, j_min+1);
            else
                c = LARGE_NUMBER;
            end
            if j_min > 1
                d = dist_fct_matrix(i_min-1, j_min-1);
            else
                d = LARGE_NUMBER;
            end
            x0 = dist_fct_matrix(i_min, j_min);
            trial_grid_points(num_trail,3) = fzero(@eikonal,x0,[],a,b,c,d,h); 
            dist_fct_matrix(i_min-1, j_min) = trial_grid_points(num_trail,3);
        end
    end
    % update the value of the upper neighbour
    if j_min < num_y_grid_lines
        if status_matrix(i_min, j_min+1) ~= 2
            if status_matrix(i_min, j_min+1) == 0
                % this point is not yet in the trail list, so put it there
                num_trail = num_trail+1;
                status_matrix(i_min, j_min+1) = 1;

                trial_grid_points(num_trail,1) = i_min;
                trial_grid_points(num_trail,2) = j_min+1;
            end
            
            % get the new dist fct value of this trail point
            if i_min < num_x_grid_lines
                a = dist_fct_matrix(i_min+1  ,j_min+1);
            else 
                a = LARGE_NUMBER;
            end
            if i_min > 1
                b = dist_fct_matrix(i_min-1  ,j_min+1);
            else
                b = LARGE_NUMBER;
            end
            if j_min < num_y_grid_lines-1
                c = dist_fct_matrix(i_min    ,j_min+2);
            else
                c = LARGE_NUMBER;
            end
            d = dist_fct_matrix(i_min    ,j_min);
            x0 = dist_fct_matrix(i_min   ,j_min);
            trial_grid_points(num_trail,3) = fzero(@eikonal,x0,[],a,b,c,d,h);
            dist_fct_matrix(i_min, j_min+1) = trial_grid_points(num_trail,3);
        end
    end
    % update the value of the lower neighbour
    if j_min > 1
        if status_matrix(i_min, j_min-1)  ~= 2
            if status_matrix(i_min, j_min-1) == 0
                % this point is not yet in the trail list, so put it there
                num_trail = num_trail+1;
                status_matrix(i_min, j_min-1) = 1;

                trial_grid_points(num_trail,1) = i_min;
                trial_grid_points(num_trail,2) = j_min-1;
            end
            
            % get the new dist fct value of this trail point
            if i_min < num_x_grid_lines
                a = dist_fct_matrix(i_min+1,j_min-1);
            else
                a = LARGE_NUMBER;
            end
            if i_min > 1
                b = dist_fct_matrix(i_min-1,j_min-1);
            else
                b = LARGE_NUMBER;
            end
            c = dist_fct_matrix(i_min,j_min);
            if j_min > 2
                d = dist_fct_matrix(i_min,j_min-2);
            else 
                d = LARGE_NUMBER;
            end
            x0 = dist_fct_matrix(i_min,j_min);
            trial_grid_points(num_trail,3) = fzero(@eikonal,x0,[],a,b,c,d,h);   
            dist_fct_matrix(i_min, j_min-1) = trial_grid_points(num_trail,3);
        end
    end
end % while there exist trail points 


function f = eikonal(x,a,b,c,d,h)
% a: right
% b: left
% c: upper
% d: lower

Dx_minus    = (x-b)/h;
Dx_plus     = (a-x)/h;
Dy_minus    = (x-d)/h;
Dy_plus     = (c-x)/h;

f = (max([Dx_minus, -Dx_plus, 0]))^2 + (max([Dy_minus, -Dy_plus, 0]))^2 - 1; 





