function [domain] = lsGenerateGridLines(domain, sc)

% sc: supersampling coefficient

index = 1;
for i_x = 1: domain.x_spacing : domain.x_size
    domain.x_grid_lines(index) = i_x;
    index = index+1;   
end

index = 1;
for i_y = 1:domain. y_spacing : domain.y_size
    domain.y_grid_lines(index) = i_y;   
    index = index + 1;
end

% generate fine grid linesindex = 1;
for i_x = 1: domain.x_spacing/sc : domain.x_size
    domain.x_grid_lines_f(index) = i_x;
    index = index+1;   
end

index = 1;
for i_y = 1:domain. y_spacing/sc : domain.y_size
    domain.y_grid_lines_f(index) = i_y;   
    index = index + 1;
end