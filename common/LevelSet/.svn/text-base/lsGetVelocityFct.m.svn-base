function F = lsGetVelocityFct(phi, phi_target, i_end, j_end, grad_x, grad_y, domain)
% LSGETVELOCITYFCT calculates the speed F at the grid points
%    
%
%
% SYNOPSIS   F = lsGetVelocityFct(phi, phi_target, i_end, j_end, grad_x, grad_y, domain)
%
%
% INPUT      phi        :
%            phi_target :
%            i_end      :
%            j_end      :
%            grad_x     :
%            grad_y     :
%            domain     :
%                          
% 
% OUTPUT     F         :
%              
%                           
% DEPENDENCES    lsGetVelocityFct uses {                                
%                                       }
%
%                lsGetVelocityFct is used by { 
%                                           }
%
% Matthias Machacek 06/24/04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculate the speed function F  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
   % Uniform velocity field
   for i=1:i_end
      for j=1:j_end
         F(i,j) = 1;
      end
   end
elseif 0
   % Flow in x-direction, y- direction, circular velocity field
   for i=1:i_end
      for j=1:j_end
         grad_l = sqrt(grad_x(i,j)^2 + grad_y(i,j)^2);
         if grad_l == 0
             grad_l = 0.000001;
         end
         
         % x-translation
         %F(i,j) =  grad_x(i,j) / grad_l;
         
         % y-translation
         % F(i,j) =  grad_y(i,j) / grad_l;
         
         % rotation
         F(i,j) =  -j/2*grad_x(i,j) / grad_l + i/2*grad_y(i,j) / grad_l;
      end
   end
elseif 0
   % Gradient driven flow
   for i=1:i_end
      for j=1:j_end
         grad_phi_x = max(A(i,j),0) + min(B(i,j),0);
         grad_phi_y = max(C(i,j),0) + min(D(i,j),0);
         grad_phi_l = sqrt(grad_phi_x^2 + grad_phi_y^2);
         
         F(i,j) =  -Dx_plus_target(i,j) * grad_x - Dy_plus_target(i,j) * grad_y;
      end
   end  
elseif 1
    % Distance driven flow
    vel_m = 5;
    if vel_m == 1
        F = phi - phi_target;
    elseif vel_m == 2
        F = sign(phi - phi_target);
    elseif vel_m == 3
        F = atan(phi - phi_target);
    elseif vel_m == 4
        F = asinh(phi - phi_target);    
    elseif vel_m == 5
        kappa = lsCurvature(phi, domain.x_spacing, domain.y_spacing, i_end, j_end);
        d_level = phi - phi_target;
        d_l = d_level >= 0;
        % protrusion
        F_prot = asinh(d_level.*d_l);
        % retraction
        F_ret = d_level.*(~d_l) .* (0.3 + kappa);
        %F_ret = d_level.*(~d_l) .* (10.3 + kappa);
        F= F_prot + F_ret;
    end
else
    % Curvature driven flow
    F = lsCurvature(phi, domain.x_spacing, domain.y_spacing, i_end, j_end);
    %F = 1 - 0.35 .* F;
    F = - F;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%