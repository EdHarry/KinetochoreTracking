function x_velocity = ls_dis_interpolate(F, x, domain, i_end, j_end)
% interpolate F at the the positions x

num_points = size(x,1);


delta_x_inv = 1/domain.x_spacing;
delta_y_inv = 1/domain.y_spacing;

for i = 1:num_points
    % find the corners belonging to x(i)
    j=1;
    while x(i,1) >= domain.x_grid_lines(j)
        j=j+1;
    end
    
    x_2 = j;
    x_1 = j-1;
    
    j=1;
    while x(i,2) >= domain.y_grid_lines(j)
        j=j+1;
    end
    
    y_2 = j;
    y_1 = j-1;  
    

    %
    %       |      |
    %-------1------2--------
    %       |      |
    %       |      |
    %-------4------3-------- ->x
    %       |      |
    %              \/y        
    
    % check for evolving edges in the x direction
    if x_1 > 1 & x_2 < i_end-1 & y_1 > 1 & y_2 < j_end-1
        Dx_minus_1 = (F(x_1  , y_1) - F(x_1-1, y_1)) * delta_x_inv;
        Dx_plus_1  = (F(x_1+1, y_1) - F(x_1  , y_1)) * delta_x_inv;
        Dx_plus_2  = (F(x_1+2, y_1) - F(x_1+1, y_1)) * delta_x_inv;
        Dx_minus_2 = Dx_plus_1;

        Dx_minus_4 = (F(x_1  , y_1+1) - F(x_1-1, y_1+1)) * delta_x_inv;
        Dx_plus_4  = (F(x_1+1, y_1+1) - F(x_1  , y_1+1)) * delta_x_inv;
        Dx_plus_3  = (F(x_1+2, y_1+1) - F(x_1+1, y_1+1)) * delta_x_inv;
        Dx_minus_3 = Dx_plus_4;

        % check for evolving edges in the y direction
        Dy_minus_1 = (F(x_1, y_1  ) - F(x_1, y_1-1)) * delta_y_inv;
        Dy_plus_1  = (F(x_1, y_1+1) - F(x_1, y_1  )) * delta_y_inv;
        Dy_plus_4  = (F(x_1, y_1+2) - F(x_1, y_1+1)) * delta_y_inv;
        Dy_minus_4 = Dy_plus_1;

        Dy_minus_2 = (F(x_1+1, y_1  ) - F(x_1+1, y_1-1)) * delta_y_inv;
        Dy_plus_2  = (F(x_1+1, y_1+1) - F(x_1+1, y_1  )) * delta_y_inv;
        Dy_plus_3  = (F(x_1+1, y_1+2) - F(x_1+1, y_1+1)) * delta_y_inv;
        Dy_minus_3 = Dy_plus_2;

        dx12 = Dx_minus_1 * Dx_plus_2;
        dx43 = Dx_minus_4 * Dx_plus_3;
        dy14 = Dy_minus_1 * Dy_plus_4;
        dy23 = Dy_minus_2 * Dy_plus_3;
    
        if (dy14 < 0) & (dy23 < 0)
            % case 1
            % upper
            x_c = [domain.x_grid_lines(x_1);domain.x_grid_lines(x_1);domain.x_grid_lines(x_2);domain.x_grid_lines(x_2)];
            y_c = [domain.y_grid_lines(y_1-1);domain.y_grid_lines(y_1);domain.y_grid_lines(y_1-1);domain.y_grid_lines(y_1)];
            i_values = [F(x_1, y_1-1);F(x_1, y_1);F(x_2, y_1-1);F(x_2, y_1)];
            i_points = [x_c,y_c];
            st = tpaps(i_points',i_values',0);
            x_velocity(i)=fnval(st,x(i,:)');
            plot(x(i,1),x(i,2),'.r');
            
        elseif (dx12 < 0) & (dx43 < 0)
            % normal pushed interpolation
            x_c = [domain.x_grid_lines(x_1);domain.x_grid_lines(x_1);domain.x_grid_lines(x_2);domain.x_grid_lines(x_2)];
            y_c = [domain.y_grid_lines(y_1);domain.y_grid_lines(y_2);domain.y_grid_lines(y_1);domain.y_grid_lines(y_2)];
            i_values = [F(x_1, y_1);F(x_1, y_2);F(x_2, y_1);F(x_2, y_2)];
            i_points = [x_c,y_c];
            st = tpaps(i_points',i_values',0);
            x_velocity(i)=1.4 .* fnval(st,x(i,:)');

            % case 2
            % left
%             x_c = [domain.x_grid_lines(x_1-1);domain.x_grid_lines(x_1-1);domain.x_grid_lines(x_1);domain.x_grid_lines(x_1)];
%             y_c = [domain.y_grid_lines(y_1);domain.y_grid_lines(y_2);domain.y_grid_lines(y_1);domain.y_grid_lines(y_2)];
%             i_values = [F(x_1-1, y_1);F(x_1-1, y_2);F(x_1, y_1);F(x_1, y_2)];
%             i_points = [x_c,y_c];
%             st = tpaps(i_points',i_values',0);
            
%             % right 
%             x_c = [domain.x_grid_lines(x_2);domain.x_grid_lines(x_2);domain.x_grid_lines(x_2+1);domain.x_grid_lines(x_2+1);domain.x_grid_lines(x_2+2);domain.x_grid_lines(x_2+2)];
%             y_c = [domain.y_grid_lines(y_1);domain.y_grid_lines(y_2);domain.y_grid_lines(y_1);domain.y_grid_lines(y_2);domain.y_grid_lines(y_1);domain.y_grid_lines(y_2)];
%             i_values = [F(x_2, y_1);F(x_2, y_2);F(x_2+1, y_1);F(x_2+1, y_2);F(x_2+2, y_1);F(x_2+2, y_2)];
%             %i_values = [F(y_1, x_2);F(y_2, x_2);F(y_1, x_2+1);F(y_2, x_2+1);F(x_1, x_2+2);F(y_2, x_2+2)];
%             i_points = [x_c,y_c];
%             st = tpaps(i_points',i_values',1);            
%             x_velocity(i)=fnval(st,x(i,:)');
%             plot(x(i,1),x(i,2),'*g');
%             [X,Y] = meshgrid(domain.x_grid_lines, domain.y_grid_lines);
%             
%             surface(X, Y, F');
%             plot(i_points(:,1), i_points(:,2),'r.');
%             %plot3(i_points(:,1), i_points(:,2),i_values, 'r.');
            
        elseif (dx12 < 0) & (dy23 < 0)
            % case 3  
            % take 1,3,4
            x_c = [domain.x_grid_lines(x_1);domain.x_grid_lines(x_2);domain.x_grid_lines(x_1)];
            y_c = [domain.y_grid_lines(y_1);domain.y_grid_lines(y_2);domain.y_grid_lines(y_2)];
            i_values = [F(x_1, y_1);F(x_2, y_2);F(x_1, y_2)];
            i_points = [x_c,y_c];
            st = tpaps(i_points',i_values',0);
            x_velocity(i)=fnval(st,x(i,:)');
            plot(x(i,1),x(i,2),'.r');
            
        elseif (dx43 < 0) & (dy14 < 0)
            % case 4
            % take 1,2,3
            x_c = [domain.x_grid_lines(x_1);domain.x_grid_lines(x_2);domain.x_grid_lines(x_2)];
            y_c = [domain.y_grid_lines(y_1);domain.y_grid_lines(y_1);domain.y_grid_lines(y_2)];
            i_values = [F(x_1, y_1);F(x_2, y_1);F(x_2, y_2)];
            i_points = [x_c,y_c];
            st = tpaps(i_points',i_values',0);
            x_velocity(i)=fnval(st,x(i,:)');
            plot(x(i,1),x(i,2),'.y');
            
        elseif (dx12 < 0) & (dy14 < 0)
            % case 5
            % take 2,3,4
            x_c = [domain.x_grid_lines(x_2);domain.x_grid_lines(x_2);domain.x_grid_lines(x_1)];
            y_c = [domain.y_grid_lines(y_1);domain.y_grid_lines(y_2);domain.y_grid_lines(y_2)];
            i_values = [F(x_2, y_1);F(x_2, y_2);F(x_1, y_2)];
            i_points = [x_c,y_c];
            st = tpaps(i_points',i_values',0);
            x_velocity(i)=fnval(st,x(i,:)');
            plot(x(i,1),x(i,2),'.c');
            
        elseif (dx43 < 0) & (dy23 < 0)
            % case 6
            % take 1,2,4
            x_c = [domain.x_grid_lines(x_1);domain.x_grid_lines(x_2);domain.x_grid_lines(x_1)];
            y_c = [domain.y_grid_lines(y_1);domain.y_grid_lines(y_1);domain.y_grid_lines(y_2)];
            i_values = [F(x_1, y_1);F(x_2, y_1);F(x_1, y_2)];
            i_points = [x_c,y_c];
            st = tpaps(i_points',i_values',0);
            x_velocity(i)=fnval(st,x(i,:)');
            plot(x(i,1),x(i,2),'.y');
            
        else
            % normal interpolation
            x_c = [domain.x_grid_lines(x_1);domain.x_grid_lines(x_1);domain.x_grid_lines(x_2);domain.x_grid_lines(x_2)];
            y_c = [domain.y_grid_lines(y_1);domain.y_grid_lines(y_2);domain.y_grid_lines(y_1);domain.y_grid_lines(y_2)];
            i_values = [F(x_1, y_1);F(x_1, y_2);F(x_2, y_1);F(x_2, y_2)];
            i_points = [x_c,y_c];
            st = tpaps(i_points',i_values',0);
            x_velocity(i)=fnval(st,x(i,:)');
            
            %velocity_fct_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, F');
            %x_velocity = fnval(velocity_fct_spline, x');
        end
    else
        
        % normal interpolation
        x_c = [domain.x_grid_lines(x_1);domain.x_grid_lines(x_1);domain.x_grid_lines(x_2);domain.x_grid_lines(x_2)];
        y_c = [domain.y_grid_lines(y_1);domain.y_grid_lines(y_2);domain.y_grid_lines(y_1);domain.y_grid_lines(y_2)];
        i_values = [F(x_1, y_1);F(x_1, y_2);F(x_2, y_1);F(x_2, y_2)];
        i_points = [x_c,y_c];
        st = tpaps(i_points',i_values',0);
        x_velocity(i)=fnval(st,x(i,:)');

    end
        % calculate weights
%         h1 = sqrt((x(i,1)-domain.x_grid_lines(x_1))^2 + (x(i,2)-domain.y_grid_lines(y_1))^2);
%         h2 = sqrt((x(i,1)-domain.x_grid_lines(x_2))^2 + (x(i,2)-domain.y_grid_lines(y_1))^2);
%         h3 = sqrt((x(i,1)-domain.x_grid_lines(x_2))^2 + (x(i,2)-domain.y_grid_lines(y_2))^2);        
%         h4 = sqrt((x(i,1)-domain.x_grid_lines(x_1))^2 + (x(i,2)-domain.y_grid_lines(y_2))^2); 
%         p = 2;
%         ht = h1^p + h2^p + h3^p + h4^p;
%         w1 = h1^p / ht;
%         w2 = h2^p / ht;
%         w3 = h3^p / ht;
%         w4 = h4^p / ht;
%         
%         %calculate nodal functions
%         Q1 = Dx_plus_1*(x(i,1)-domain.x_grid_lines(x_1)) + Dy_plus_1*(x(i,2)-domain.y_grid_lines(y_1)) +  F(x_1, y_1);
%         Q2 = Dx_plus_2*(x(i,1)-domain.x_grid_lines(x_2)) + Dy_plus_2*(x(i,2)-domain.y_grid_lines(y_1)) +  F(x_2, y_1);
%         Q3 = Dx_plus_3*(x(i,1)-domain.x_grid_lines(x_2)) + Dy_plus_3*(x(i,2)-domain.y_grid_lines(y_2)) +  F(x_2, y_2);
%         Q4 = Dx_plus_4*(x(i,1)-domain.x_grid_lines(x_1)) + Dy_plus_4*(x(i,2)-domain.y_grid_lines(y_2)) +  F(x_1, y_2);
%         
%         x_velocity(i) = w1*Q1 + w2*Q2 + w3*Q3 + w4*Q4;
        
%     end
    
end
