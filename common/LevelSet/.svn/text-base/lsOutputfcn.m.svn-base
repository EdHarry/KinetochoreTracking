%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = lsOutputfcn(t,y, flag, phi_t1, i_end, j_end, delta_x, delta_y, domain)
f = figure(gcf);

if ~isempty(y)
    if 1
        phi_vec = y(1:i_end*j_end, end);
        phi = reshape(phi_vec , i_end, j_end);
        x_vec = y(i_end*j_end+1 : end, end);
        x = reshape(x_vec, length(x_vec)/2,2);
        % Extract the zero level
        phi_zero = lsGetZeroLevel(phi, domain);

        plot(phi_zero(1,:), phi_zero(2,:),'k');  
        
        % plot the tracking points 
        % plot(x_vec(1,:),x_vec(2,:),'.r');
    else
        phi = reshape(y, i_end, j_end);
        % Extract the zero level
        phi_zero = lsGetZeroLevel(phi(:,:,end), domain);
        plot(phi_zero(1,:), phi_zero(2,:),'k');       
    end
end

status = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%