function [phi_next, F, delta_plus, delta_minus, dt] = lsSolveConvection(phi, phi_p, F_p, delta_plus_p, delta_minus_p, delta_t, delta_x, delta_y, i_end, j_end, phi_target, domain)
% LSSOLVECONVECTION solves the equation d phi/ dt + F|grad phi| = 0 
%    
%           Second order spatial
%           and first order temporal scheme 
%           according to J.A. Sethian: Level Set Methods
%           and Fast Marching Methods p.: 66
%
%
% SYNOPSIS   [phi_next, velocity_fct, delta_t_optimal] = lsSolveConvection(phi, delta_t, delta_x, delta_y, i_end, j_end, phi_target, domain)
%
%
% INPUT      phi        :   function value at current time step t
%            phi_p      :   function value at previous time step t-1
%            F_p        :   velocity at previous time step
%            delta_t    :
%            delta_x    :   
%            delta_y    :
%            i_end      :
%            j_end      :
%            phi_target :
%            domain     :
%                          
% 
% OUTPUT     phi_next       : function value at next time step t+1
%            velocity_fct   : used velocity fct
%            dt             : used time step 
%              
%                           
% DEPENDENCES    lsSolveConvection uses {                                
%                                       }
%
%                lsSolveConvection is used by { 
%                                           }
%
% Matthias Machacek 06/22/04


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%% Difference operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(phi, delta_x, delta_y, i_end, j_end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% figure
% quiver(domain.x_grid_lines, domain.y_grid_lines, grad_x, grad_y);
% hold on
% surface(domain.x_grid_lines,domain.y_grid_lines,phi);
% contour(domain.x_grid_lines,domain.y_grid_lines,phi, 40);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Difference operators for the target function %%%%%%
for i=1:i_end
   for j=1:j_end
      if i > i_end-1  
         Dx_plus_target(i,j) = (phi_target(i  ,j)-phi_target(i-1,j)) / delta_x;  
      else
         Dx_plus_target(i,j) = (phi_target(i+1,j)-phi_target(i  ,j)) / delta_x;
      end
      
      if j > j_end-1
         Dy_plus_target(i,j) = (phi_target(i,j  )-phi_target(i,j-1)) / delta_y;
      else
         Dy_plus_target(i,j) = (phi_target(i,j+1)-phi_target(i,j  )) / delta_y;
      end
   end
end

% get the gradient at the zero level set position of the target function
for i=1:i_end
   for j=1:j_end
      grad_phi_target(i,j) = sqrt(Dx_plus_target(i,j)^2+Dy_plus_target(i,j)^2);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculate the speed function F  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = lsGetVelocityFct(phi, phi_target, i_end, j_end, grad_x, grad_y, domain);


% get the maximal time step based on the CLF number 
delta_t_optimal = 0.9 * delta_x / max(max(F));
%dt = delta_t_optimal;
dt = delta_t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update the level-sets   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% methods Euler, Adams, Heun

method = 'Heun';  

if strcmp(method,'Euler')
    % First order acccurate: Euler method
    for i=1:i_end
        for j=1:j_end
            phi_next(i,j) = phi(i,j) - dt*(max(F(i,j), 0)*delta_plus(i,j) + min(F(i,j), 0) * delta_minus(i,j));
        end
    end
elseif strcmp(method,'Adams')
    % Adams_Bashforth method: second order
    for i=1:i_end
        for j=1:j_end
            a = max(F(i,j), 0)*delta_plus(i,j) + min(F(i,j), 0) * delta_minus(i,j);
            b = max(F_p(i,j), 0)*delta_plus_p(i,j) + min(F_p(i,j), 0) * delta_minus_p(i,j);

            phi_next(i,j) = phi(i,j) - dt*(3/2  * a - 0.5 * b);
        end
    end
elseif strcmp(method,'Heun')
    % Heun method: second order. A predictor-corrector schema
    % Predictor step
    for i=1:i_end
        for j=1:j_end
            phi_p(i,j) = phi(i,j) - dt*(max(F(i,j), 0)*delta_plus(i,j) + min(F(i,j), 0) * delta_minus(i,j));
        end
    end
   
   [delta_plus_p, delta_minus_p, grad_x_p, grad_y_p] = lsGradient2o(phi_p,...
                                           delta_x, delta_y, i_end, j_end);
   
   F_p = lsGetVelocityFct(phi, phi_target, i_end, j_end, grad_x_p, grad_y_p, domain);
   
   % Corrector step
   for i=1:i_end
      for j=1:j_end
         phi_next(i,j) = phi(i,j) -...
             dt/2*(max(F(i,j),   0)*delta_plus(i,j)   + min(F(i,j),   0) * delta_minus(i,j)+...
                   max(F_p(i,j), 0)*delta_plus_p(i,j) + min(F_p(i,j), 0) * delta_minus_p(i,j));
      end
   end   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%