function track_points = lsIntegrateVelocity(phi, F, t_steps,... 
                     grid_coordinates, delta_t, delta_x, delta_y, i_end, j_end, domain)
% LSINTEGRATEVELOCITY integrates velocity to get time dependent position
%    
%
%
% SYNOPSIS   track_points = lsIntegrateVelocity(dist_matrix, velocity_fct_matrix, grid_coordinates, delta_t, delta_x, delta_y, i_end, j_end, domain)
%
%
% INPUT      dist_matrix            :
%            velocity_fct_matrix    :
%            grid_coordinates       :
%            delta_t                :   
%            delta_x                :
%            delta_y                :
%            i_end                  :
%            j_end                  :
%            domain                 :
%                          
% 
% OUTPUT     track_points           :
%              
%                           
% DEPENDENCES    lsIntegrateVelocity uses {                                
%                                       }
%
%                lsIntegrateVelocity is used by { 
%                                           }
%
% Matthias Machacek 06/24/04

contr = 0;

track_points(:,:,1) = lsGetZeroLevel(phi(:,:,1), domain);

% number of time points 
num_time_steps = size(phi,3)-1;


for i = 1:num_time_steps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Velocity interpolation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find a B-spline interpolation of the velocity field
    velocity_fct_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, F(:,:,i)');
    
    % Get velocity at these points 
    track_points_velocity = fnval(velocity_fct_spline, track_points(:,:,i));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Gradient field interpolation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the gradient at these points (grad phi)
    if 1
        [grad_x, grad_y] = lsGradient(phi(:,:,i), 2, 0, delta_x, delta_y, i_end, j_end);
        
        % Find a B-spline interpolation of the gradient field
        grad_x_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_x');
        grad_y_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_y');

        % Get the gradient at the track points
        track_points_grad_x = fnval(grad_x_spline, track_points(:,:,i));
        track_points_grad_y = fnval(grad_y_spline, track_points(:,:,i));

        grad = sqrt(track_points_grad_x.^2 + track_points_grad_y.^2);
        
        track_points_grad_x_u = track_points_grad_x ./ grad;
        track_points_grad_y_u = track_points_grad_y ./ grad;      
        
    elseif 0
        [grad_x_l, grad_y_l] = lsGradient(phi(:,:,i), 1,  -1, delta_x, delta_y, i_end, j_end);
        [grad_x_r, grad_y_r] = lsGradient(phi(:,:,i), 1,  1, delta_x, delta_y, i_end, j_end);
        
        % Find a B-spline interpolation of the gradient field
        grad_x_l_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_x_l');
        grad_y_l_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_y_l');
        grad_x_r_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_x_r');
        grad_y_r_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_y_r');
        
        % Get the gradient at the track points
        track_points_grad_x_l = fnval(grad_x_l_spline, track_points(:,:,i));
        track_points_grad_y_l = fnval(grad_y_l_spline, track_points(:,:,i));
        track_points_grad_x_r = fnval(grad_x_r_spline, track_points(:,:,i));
        track_points_grad_y_r = fnval(grad_y_r_spline, track_points(:,:,i));
        
        % Summation
        grad_n_rr = sqrt(track_points_grad_x_r.^2 + track_points_grad_y_r.^2); 
        grad_n_rl = sqrt(track_points_grad_x_r.^2 + track_points_grad_y_l.^2);        
        grad_n_lr = sqrt(track_points_grad_x_l.^2 + track_points_grad_y_r.^2);     
        grad_n_ll = sqrt(track_points_grad_x_l.^2 + track_points_grad_y_l.^2);    
        
        grad_x_s = track_points_grad_x_r ./ grad_n_rr +...
                   track_points_grad_x_r ./ grad_n_rl +...
                   track_points_grad_x_l ./ grad_n_lr +...
                   track_points_grad_x_l ./ grad_n_ll; 
               
        grad_y_s = track_points_grad_y_r ./ grad_n_rr +...
                   track_points_grad_y_l ./ grad_n_rl +...
                   track_points_grad_y_r ./ grad_n_lr +...
                   track_points_grad_y_l ./ grad_n_ll; 
               

    
        grad = sqrt(grad_x_s.^2 + grad_y_s.^2);
    
        track_points_grad_x_u = grad_x_s ./ grad;
        track_points_grad_y_u = grad_y_s ./ grad;
    else
        % Use the matlab build-in gradient function
        [grad_x, grad_y] = gradient(phi(:,:,i), delta_x, delta_y);

        % Find a B-spline interpolation of the gradient field
        grad_x_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_x');
        grad_y_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_y');

        % Get the gradient at the track points
        track_points_grad_x = fnval(grad_x_spline, track_points(:,:,i));
        track_points_grad_y = fnval(grad_y_spline, track_points(:,:,i));

        grad = sqrt(track_points_grad_x.^2 + track_points_grad_y.^2);

        track_points_grad_x_u =  track_points_grad_x./ grad;
        track_points_grad_y_u =  track_points_grad_y./ grad;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    if 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Integrate velocity  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        track_points(1,:,i+1) = track_points(1,:,i) +  track_points_grad_x_u .* track_points_velocity * delta_t(i);
        track_points(2,:,i+1) = track_points(2,:,i) +  track_points_grad_y_u .* track_points_velocity * delta_t(i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    elseif 1
        % Heun method: second order. A predictor-corrector schema
        
        %%%% Predictor step  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        track_points_p(1,:) = track_points(1,:,i) +  track_points_grad_x_u .* track_points_velocity * delta_t(i);
        track_points_p(2,:) = track_points(2,:,i) +  track_points_grad_y_u .* track_points_velocity * delta_t(i);

        %%%% Corrector step  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get the particle velocity at the predictor place  
        track_points_velocity_p = fnval(velocity_fct_spline, track_points_p);
    
        % Get the gradient at the track points
        track_points_grad_x_p = fnval(grad_x_spline, track_points_p);
        track_points_grad_y_p = fnval(grad_y_spline, track_points_p);

        grad_p = sqrt(track_points_grad_x_p.^2 + track_points_grad_y_p.^2);

        track_points_grad_x_u_p =  track_points_grad_x_p./ grad_p;
        track_points_grad_y_u_p =  track_points_grad_y_p./ grad_p;
        
        
        track_points(1,:,i+1) = track_points(1,:,i) + delta_t(i)/2*...
            (track_points_grad_x_u   .* track_points_velocity +...
             track_points_grad_x_u_p .* track_points_velocity);
        
        track_points(2,:,i+1) = track_points(2,:,i) + delta_t(i)/2*... 
            (track_points_grad_y_u   .* track_points_velocity +...
             track_points_grad_y_u_p .* track_points_velocity_p);

    else
        % Adams_Bashforth method: second order
        
        a_x = track_points_grad_x_u   .* track_points_velocity;
        a_y = track_points_grad_y_u   .* track_points_velocity; 
        
        if i == 1
            b_x = a_x;
            b_y = a_y;           
        else
            b_x = a_x_o;
            b_y = a_y_o;
        end
        
        track_points(1,:,i+1) = track_points(1,:,i) + delta_t(i) * (3/2  * a_x - 0.5 * b_x);
        track_points(2,:,i+1) = track_points(2,:,i) + delta_t(i) * (3/2  * a_y - 0.5 * b_y);
        
        a_x_o = a_x;
        a_y_o = a_y;
    end
end







% not used control plotting
% if contr
%     figure
%     fnplt(velocity_fct_spline), axis equal
%     hold on
%     plot3(track_points(1,:,i), track_points(2,:,i), track_points_velocity, 'ro');
% end
% if contr
%     [m_grad_x, m_grad_y] = gradient(dist_matrix(:,:,i), delta_x, delta_y);
% 
%     figure
%     quiver(domain.x_grid_lines, domain.y_grid_lines, grad_x, grad_y,0);
%     hold on
%     contour(domain.x_grid_lines, domain.y_grid_lines, dist_matrix(:,:,i));
%     %quiver(domain.x_grid_lines, domain.y_grid_lines, m_grad_x, m_grad_y);
%     plot(track_points(1,:,i), track_points(2,:,i),'r');
%     %axis equal
%     %quiver(domain.y_grid_lines, domain.x_grid_lines, grad_y, grad_x);
%     quiver(track_points(1,:,i), track_points(2,:,i), track_points_grad_x, track_points_grad_y,0,'r');
%     axis equal
%     xlabel('x');
%     ylabel('y');
% 
%     for ii=1:size(grad_x,1)
%         for jj=1:size(grad_x,2)
%             grad_x_vec((ii-1)*size(grad_x,2)+jj) = grad_x(ii,jj);
%             grad_y_vec((ii-1)*size(grad_y,2)+jj) = grad_y(ii,jj);
%             x_cord((ii-1)*size(grad_x,2)+jj) = domain.x_grid_lines(ii);
%             y_cord((ii-1)*size(grad_x,2)+jj) = domain.y_grid_lines(jj);
%         end
%     end
% 
%     figure
%     fnplt(grad_x_spline)
%     hold on
%     plot3(track_points(1,:), track_points(2,:), track_points_grad_x,'o');
%     plot3(x_cord, y_cord, grad_x_vec, 'x');
%     xlabel('x');
%     ylabel('y');
% 
%     figure
%     fnplt(grad_y_spline)
%     hold on
%     plot3(track_points(1,:), track_points(2,:), track_points_grad_y,'o');
%     plot3(x_cord, y_cord, grad_y_vec, 'x');
%     xlabel('x');
%     ylabel('y');
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 





% options = [];
% x_0 = reshape(x_0, prod(size(x_0)), 1);
% global tt;
% tt=1;
% [t_steps, Y, TE,YE,IE] = ode45(@dx_fct, t_steps, x_0,  options,...
%     phi, F, i_end, j_end, delta_x, delta_y, domain);
% 
% 
% track_points = reshape(Y, length(Y)/2, 2);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function dx_vec = dx_fct(t, y, phi, F, i_end, j_end, delta_x, delta_y, domain)
% global tt;
% x = reshape(y, length(y)/2, 2);
% 
% time_step = length(t);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Velocity interpolation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % find a B-spline interpolation of the velocity field
% velocity_fct_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, F(:,:,tt)');
% 
% % Get velocity at these points
% x_velocity = fnval(velocity_fct_spline, x');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Use the matlab build-in gradient function
% [grad_x, grad_y] = gradient(phi(:,:,tt), delta_x, delta_y);
% 
% % Find a B-spline interpolation of the gradient field
% grad_x_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_x');
% grad_y_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_y');
% 
% % Get the gradient at the track points
% track_points_grad_x = fnval(grad_x_spline, x');
% track_points_grad_y = fnval(grad_y_spline, x');
% 
% grad = sqrt(track_points_grad_x.^2 + track_points_grad_y.^2);
% 
% track_points_grad_x_u =  track_points_grad_x./ grad;
% track_points_grad_y_u =  track_points_grad_y./ grad;
% tt=tt+1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
% dx_x = x_velocity .* track_points_grad_x./ grad;
% dx_y = x_velocity .* track_points_grad_y./ grad;
% 
% dx_vec = cat(1,dx_x',dx_y');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




