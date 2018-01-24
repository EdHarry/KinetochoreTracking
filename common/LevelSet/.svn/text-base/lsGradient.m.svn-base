function [grad_x, grad_y]  = lsGradient(phi, order, support, delta_x, delta_y, i_end, j_end)
% LSGRADIENT calculates the gradient of order "order"
%    
%
% SYNOPSIS    [grad_x, grad_y]  = lsGradient(phi, order, support, delta_x, delta_y, i_end, j_end)
%
% INPUT         phi     :
%               order   :                                                       
%               support :    
%               delta_x : 
%               delta_y :
%               i_end   :
%               j_end   :
%                         
% 
% OUTPUT        grad_x  :
%               grad_y  :
%        
%                           
% DEPENDENCES   lsIntersectSplineLine uses {                            
%                                       }
%
%               lsIntersectSplineLine is used by { 
%                                           }
%
% Matthias Machacek 07/12/04


%%%%%%%%%%%%%%% Difference operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if order == 1
    if support == -1
        %%%%%%%%%%%%%%% First order operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dy_minus =  zeros(i_end, j_end);
        Dx_minus =  zeros(i_end, j_end);

        Dx_minus(:,2:j_end) = diff(phi,order,2);
        Dx_minus(:,1) = Dx_minus(:,2);
        Dx_minus = Dx_minus ./ delta_x;

        Dy_minus(2:i_end,:) = diff(phi,order,1);
        Dy_minus(1,:) = Dy_minus(1,:);
        Dy_minus = Dy_minus ./ delta_y;

        grad_x = Dx_minus;
        grad_y = Dy_minus;
    elseif support == 1
        Dx_plus  =  zeros(i_end, j_end);
        Dy_plus  =  zeros(i_end, j_end);

        Dx_plus = diff(phi,order,2);
        Dx_plus(:,j_end) = Dx_plus(:,j_end-1);
        Dx_plus = Dx_plus ./ delta_x;

        Dy_plus = diff(phi,order,1);
        Dy_plus(i_end,:) = Dy_plus(i_end-1,:);
        Dy_plus = Dy_plus ./ delta_y;

        grad_x = Dx_plus;
        grad_y = Dy_plus;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif order == 2
    Dx_central =  zeros(i_end, j_end);
    Dy_central =  zeros(i_end, j_end);
    
    delta_x_2  = 2  * delta_x;
    delta_y_2  = 2  * delta_y;
    
    %%%%%%%%%%%%%%% Second order central operators %%%%%%%%%%%%%%%%%%%%
    for i=1:i_end
        for j=1:j_end
            if i < 2
                Dy_central(i,j) =  (phi(i+2,j  )-  phi(i  ,j  ))  / delta_y_2;
            elseif i > i_end-1
                Dy_central(i,j) =  (phi(i  ,j  )-  phi(i-2,j  ))  / delta_y_2;
            else
                Dy_central(i,j) =  (phi(i+1,j  )-  phi(i-1,j  ))  / delta_y_2;
            end

            if j < 2
                Dx_central(i,j) =  (phi(i  ,j+2)-  phi(i  ,j  ))  / delta_x_2;
            elseif j > j_end-1
                Dx_central(i,j) =  (phi(i  ,j  )-  phi(i  ,j-2))  / delta_x_2;
            else
                Dx_central(i,j) =  (phi(i  ,j+1)-  phi(i  ,j-1))  / delta_x_2;
            end
        end
    end

    grad_x = Dx_central;
    grad_y = Dy_central;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif order == 4
    
    %%%%%%%%%%%%%%% Fourth order central operators %%%%%%%%%%%%%%%%%%%%
    Dx_central =  zeros(i_end, j_end);
    Dy_central =  zeros(i_end, j_end);
    
    delta_x_2  = 2  * delta_x;
    delta_y_2  = 2  * delta_y;
    delta_x_12 = 12 * delta_x;
    delta_y_12 = 12 * delta_y;
    
    for i=1:i_end
        for j=1:j_end
            if i == 2 || i == i_end-1
                % use second order central
                Dy_central(i,j) = (phi(i+1,j)-  phi(i-1,j))  / delta_y_2;
            elseif i == 1
                % use first order right
                Dy_central(i,j) = (phi(i+1,j) -  phi(i,j))  / delta_y;
            elseif i == i_end
                % use first order left
                Dy_central(i,j) = (phi(i,j) -  phi(i-1,j))  / delta_y;
            else
                Dy_central(i,j) = (-phi(i+2,j)+8*phi(i+1,j)-8*phi(i-1,j)+phi(i-2,j))  / delta_y_12;
            end

            
            if j == 2 || j == j_end-1
                % use second order central
                Dx_central(i,j) = (phi(i,j+1) - phi(i,j-1)) / delta_x_2;
            elseif j == 1
                % use first order right
                Dx_central(i,j) = (phi(i,j+1) - phi(i,j)) / delta_x;
            elseif j == j_end
                % use first order left
                Dx_central(i,j) = (phi(i,j) -  phi(i,j-1)) / delta_x;
            else
                Dx_central(i,j) = (-phi(i,j+2)+8*phi(i,j+1)-8*phi(i,j-1)+phi(i,j-2))  / delta_x_12;
            end
        end
    end

    grad_x = Dx_central;
    grad_y = Dy_central;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
elseif  order == 10
    

    Dy_minus  =  zeros(i_end, j_end);
    Dx_minus  =  zeros(i_end, j_end);

    Dx_plus   =  zeros(i_end, j_end);
    Dy_plus   =  zeros(i_end, j_end);



    Dx_minus(:,2:j_end) = diff(phi,1,2);
    Dx_minus(:,1) = Dx_minus(:,2);
    Dx_minus = Dx_minus ./ delta_x;

    Dy_minus(2:i_end,:) = diff(phi,1,1);
    Dy_minus(1,:) = Dy_minus(1,:);
    Dy_minus = Dy_minus ./ delta_y;

    Dx_plus = diff(phi,1,2);
    Dx_plus(:,j_end) = Dx_plus(:,j_end-1);
    Dx_plus = Dx_plus ./ delta_x;

    Dy_plus = diff(phi,1,1);
    Dy_plus(i_end,:) = Dy_plus(i_end-1,:);
    Dy_plus = Dy_plus ./ delta_y;

    
    n_x = Dx_plus ./ sqrt(Dx_plus.^2 + Dy_plus.^2) +...
          Dx_minus ./ sqrt(Dx_minus.^2 + Dy_plus.^2) +...
          Dx_plus ./ sqrt(Dx_plus.^2 + Dy_minus.^2) +...
          Dx_minus ./ sqrt(Dx_minus.^2 + Dy_minus.^2);


    n_y = Dy_plus ./ sqrt(Dx_plus.^2 + Dy_plus.^2) +...
          Dy_minus ./ sqrt(Dx_minus.^2 + Dy_plus.^2) +...
          Dy_plus ./ sqrt(Dx_plus.^2 + Dy_minus.^2) +...
          Dy_minus ./ sqrt(Dx_minus.^2 + Dy_minus.^2);
    
          
    grad_x = n_x ./ sqrt(n_x.^2 + n_y.^2);
    grad_y = n_y ./ sqrt(n_x.^2 + n_y.^2);   
end