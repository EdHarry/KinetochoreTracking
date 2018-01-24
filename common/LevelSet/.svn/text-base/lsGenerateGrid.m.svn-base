function [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain)

[x_grid, y_grid] = meshgrid(1 : domain.x_spacing : domain.x_size,  1 : domain.y_spacing : domain.y_size);
grid_coordinates(:, 1) = reshape(x_grid, prod(size(x_grid)),1);
grid_coordinates(:, 2) = reshape(y_grid, prod(size(y_grid)),1);
