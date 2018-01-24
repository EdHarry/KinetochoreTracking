function [phi_xx, phi_yy] = lsSecondDifference(phi, delta_x, delta_y, i_end, j_end)

phi_xx =  zeros(i_end, j_end);
phi_yy =  zeros(i_end, j_end);

delta_x_2  = delta_x ^2;
delta_y_2  = delta_y ^2;

%%%%%%%%%%%%%%% Second order central operators %%%%%%%%%%%%%%%%%%%%
for i=1:i_end
    for j=1:j_end
        if i == 1
            phi_yy(i,j) =  (phi(i+2,j) - 2*phi(i+1,j) + phi(i,j))  / delta_y_2;
        elseif i == i_end
            phi_yy(i,j) =  (phi(i,j) - 2*phi(i-1,j) + phi(i-2,j))  / delta_y_2;
        else
            phi_yy(i,j) =  (phi(i+1,j) - 2*phi(i,j) + phi(i-1,j))  / delta_y_2;
        end

        if j == 1
            phi_xx(i,j) =  (phi(i,j+2) - 2 * phi(i,j+1) + phi(i,j))  / delta_x_2;
        elseif j == j_end
            phi_xx(i,j) =  (phi(i,j) - 2 * phi(i,j-1) + phi(i,j-2))  / delta_x_2;
        else
            phi_xx(i,j) =  (phi(i,j+1) - 2 * phi(i,j) + phi(i,j-1))  / delta_x_2;
        end
    end
end
