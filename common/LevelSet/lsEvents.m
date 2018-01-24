%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value,isterminal,direction] = events(t,y, phi_t1, i_end, j_end, delta_x, delta_y, domain)
% value(i) is the value of the function. 
% isterminal(i) = 1 if the integration is to terminate at a zero of this event 
%     function and 0 otherwise.
% direction(i) = 0 if all zeros are to be computed (the default), +1 if only 
%        the zeros where the event function increases, and -1 if only the 
%        zeros where the event function decreases.   

global residual_last;
global residual;
global residual_i;

if 1
    phi_vec = y(1:i_end*j_end);
    phi_t = reshape(phi_vec , i_end, j_end);
    x_vec = y(i_end*j_end+1 : end);
    x = reshape(x_vec, length(x_vec)/2,2);
else
    phi_t = reshape(y, i_end, j_end);
end

res = norm(phi_t - phi_t1, 'fro');
residual_d = residual_last - res;

if 0
    if residual_d < 2
        value = 0;
    else
        value = residual_d;
    end
else
    if res < 5 | residual_d < 0
        value = -1;
    else
        value = 1;
    end
end

residual_last = res;
residual(residual_i) = res;
residual_i = residual_i+1;
isterminal = 1;
direction = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%