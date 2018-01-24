function phi_zero = lsGetZeroLevel(phi, domain, FINE)
% LSGETZEROLEVEL finds the zero level countour in a 2D matrix
% 
%
% SYNOPSIS      phi_zero = lsGetZeroLevel(phi, domain)  
%
% INPUT         phi     : 2D matrix
%               domain  : structure with x and y grid line coordinates          
% 
% OUTPUT        phi_zero: zero level contour, 2D vector (x,y)    
%                     
%                           
% DEPENDENCES    lsLineMatching uses { contourc
%                                    } 
%
% Matthias Machacek 6/22/04

% this function gets the level at the grid lines ->
% x-levels, y-levels.
% it finds them by linear interpolation of the level set 
% matrix
% interpolate Level set
if ~exist('FINE','var')
    FINE = 0;
end



if FINE
    % interpolate on finer grid
    phi_fine = interp2(domain.x_grid_lines, domain.y_grid_lines', phi,...
                    domain.x_grid_lines_f, domain.y_grid_lines_f','cubic');
end
% extract zero level set
if FINE
    phi_zero = contourc(domain.x_grid_lines_f, domain.y_grid_lines_f, phi_fine,[0 0]);
else
    phi_zero = contourc(domain.x_grid_lines, domain.y_grid_lines, phi,[0 0]);
end

if ~isempty(phi_zero)
    phi_zero(:,1)=[];
end



