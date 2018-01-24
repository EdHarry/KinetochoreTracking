function kappa = lsCurvature(phi, delta_x, delta_y, i_end, j_end)


order = 2;
[phi_x, phi_y]   = lsGradient(phi, order, 0, delta_x, delta_y, i_end, j_end);
[phi_xx, phi_yy] = lsSecondDifference(phi, delta_x, delta_y, i_end, j_end);

[phi_xx_d, phi_xy] = lsGradient(phi_x, order, 0, delta_x, delta_y, i_end, j_end);

phi_x2 = phi_x .^ 2;
phi_y2 = phi_y .^ 2;


kappa = (phi_xx .* phi_y2 - 2 .* phi_y .* phi_x .* phi_xy + phi_yy .* phi_x2) ./...
        (phi_x2 + phi_y2) .^ (3/2);
    
%kappa = (phi_xx .* phi_y2 - 2 .* phi_y .* phi_x + phi_yy .* phi_x2) ./  (phi_x2 + phi_y2) .^ (3/2);