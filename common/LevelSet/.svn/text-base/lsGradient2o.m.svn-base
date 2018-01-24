function [delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(phi, delta_x, delta_y, i_end, j_end)
% LSFIRSTSECONDDIFFERENCES calculates gradient wit hsecond order accuracy
%    
%
%
% SYNOPSIS   [delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(phi, delta_x, delta_y, i_end, j_end)
%
%
% INPUT      phi        : phi=f(x,y) function values on a grid
%            delta_x    : x-direction grid spacing
%            delta_y    : y-direction grid spacing
%            i_end      : number of x grid points
%            j_end      : number of y grid points 
%                          
% 
% OUTPUT     delta_plus     :  right side absolute value of the gradient
%            delta_minus    :  left side absolute value of the gradient
%            grad_x         :  x-comp. of the gradient
%            grad_y         :  y-comp. of the gradient
%                           
% DEPENDENCES     lsGradient2o uses {                                
%                                       }
%
%                 lsGradient2o is used by { 
%                                           }
%
% Matthias Machacek 06/22/04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%% Difference operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dy_minus          =  zeros(i_end, j_end);
Dy_plus           =  zeros(i_end, j_end);
Dy_minus_y_minus  =  zeros(i_end, j_end);
Dy_plus_y_plus    =  zeros(i_end, j_end);
Dy_plus_y_minus   =  zeros(i_end, j_end);

Dx_minus          =  zeros(i_end, j_end);
Dx_plus           =  zeros(i_end, j_end);
Dx_minus_x_minus  =  zeros(i_end, j_end);
Dx_plus_x_plus    =  zeros(i_end, j_end);
Dx_plus_x_minus   =  zeros(i_end, j_end);

delta_x_inv = 1/delta_x;
delta_y_inv = 1/delta_y;
delta_x_sq_inv = 1/delta_x^2;
delta_y_sq_inv = 1/delta_y^2;

for i=1:i_end
    if i < 3
        % assume that the value at i=1 is the same as at i=2
        Dy_minus(i,:)          =  (phi(i+1,:)-  phi(i  ,:)           )  * delta_y_inv;
        Dy_plus(i,:)           =  (phi(i+1,:)-  phi(i  ,:)           )  * delta_y_inv;
        Dy_minus_y_minus(i,:)  =  (phi(i+2,:)-2*phi(i+1,:)+phi(i,:)  )  * delta_y_sq_inv;
        Dy_plus_y_plus(i,:)    =  (phi(i+2,:)-2*phi(i+1,:)+phi(i,:)  )  * delta_y_sq_inv;
        Dy_plus_y_minus(i,:)   =  (phi(i+2,:)-2*phi(i+1,:)+phi(i,:)  )  * delta_y_sq_inv;
    elseif i > i_end-3
        Dy_minus(i,:)          =  (phi(i  ,:)-  phi(i-1,:)             )  * delta_y_inv;
        Dy_plus(i,:)           =  (phi(i  ,:)-  phi(i-1,:)             )  * delta_y_inv;
        Dy_minus_y_minus(i,:)  =  (phi(i  ,:)-2*phi(i-1,:)+phi(i-2,:)  )  * delta_y_sq_inv;
        Dy_plus_y_plus(i,:)    =  (phi(i  ,:)-2*phi(i-1,:)+phi(i-2,:)  )  * delta_y_sq_inv;
        Dy_plus_y_minus(i,:)   =  (phi(i  ,:)-2*phi(i-1,:)+phi(i-2,:)  )  * delta_y_sq_inv;
    else
        Dy_minus(i,:)          =  (phi(i  ,:)-  phi(i-1,:)             )  * delta_y_inv;
        Dy_plus(i,:)           =  (phi(i+1,:)-  phi(i  ,:)             )  * delta_y_inv;
        Dy_minus_y_minus(i,:)  =  (phi(i  ,:)-2*phi(i-1,:)+phi(i-2,:)  )  * delta_y_sq_inv;
        Dy_plus_y_plus(i,:)    =  (phi(i+2,:)-2*phi(i+1,:)+phi(i  ,:)  )  * delta_y_sq_inv;
        Dy_plus_y_minus(i,:)   =  (phi(i+1,:)-2*phi(i  ,:)+phi(i-1,:)  )  * delta_y_sq_inv;
    end
end

for j=1:j_end
    if j < 3
        Dx_minus(:,j)          =  (phi(:,j+1)-  phi(:,j  )             )  * delta_x_inv;
        Dx_plus(:,j)           =  (phi(:,j+1)-  phi(:,j  )             )  * delta_x_inv;
        Dx_minus_x_minus(:,j)  =  (phi(:,j+2)-2*phi(:,j+1)+phi(:,j  )  )  * delta_x_sq_inv;
        Dx_plus_x_plus(:,j)    =  (phi(:,j+2)-2*phi(:,j+1)+phi(:,j  )  )  * delta_x_sq_inv;
        Dx_plus_x_minus(:,j)   =  (phi(:,j+2)-2*phi(:,j+1)+phi(:,j  )  )  * delta_x_sq_inv;
    elseif j > j_end-3
        Dx_minus(:,j)          =  (phi(:,j  )-  phi(:,j-1)             )  * delta_x_inv;
        Dx_plus(:,j)           =  (phi(:,j  )-  phi(:,j-1)             )  * delta_x_inv;
        Dx_minus_x_minus(:,j)  =  (phi(:,j  )-2*phi(:,j-1)+phi(:,j-2)  )  * delta_x_sq_inv;
        Dx_plus_x_plus(:,j)    =  (phi(:,j  )-2*phi(:,j-1)+phi(:,j-2)  )  * delta_x_sq_inv;
        Dx_plus_x_minus(:,j)   =  (phi(:,j  )-2*phi(:,j-1)+phi(:,j-2)  )  * delta_x_sq_inv;
    else
        Dx_minus(:,j)          =  (phi(:,j  )-  phi(:,j-1)             )  * delta_x_inv;
        Dx_plus(:,j)           =  (phi(:,j+1)-  phi(:,j  )             )  * delta_x_inv;
        Dx_minus_x_minus(:,j)  =  (phi(:,j  )-2*phi(:,j-1)+phi(:,j-2)  )  * delta_x_sq_inv;
        Dx_plus_x_plus(:,j)    =  (phi(:,j+2)-2*phi(:,j+1)+phi(:,j  )  )  * delta_x_sq_inv;
        Dx_plus_x_minus(:,j)   =  (phi(:,j+1)-2*phi(:,j  )+phi(:,j-1)  )  * delta_x_sq_inv;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

grad_x =        zeros(i_end, j_end);
grad_y =        zeros(i_end, j_end);
delta_plus  =   zeros(i_end, j_end);
delta_minus =   zeros(i_end, j_end);

delta_x_half = delta_x/2;
delta_y_half = delta_y/2;


A1  = (Dx_minus_x_minus .* Dx_plus_x_minus) >= 0;
A2  = abs(Dx_minus_x_minus) <= abs(Dx_plus_x_minus);
A = Dx_minus + delta_x_half * (A1.*A2 .* Dx_minus_x_minus + (A1&~A2) .* Dx_plus_x_minus);

B1 = (Dx_plus_x_plus .* Dx_plus_x_minus) >= 0;
B2 = abs(Dx_plus_x_plus) <= abs(Dx_plus_x_minus);
B = Dx_plus - delta_x_half * (B1.*B2 .* Dx_plus_x_plus   + (B1&~B2) .* Dx_plus_x_minus);

C1 = (Dy_minus_y_minus .* Dy_plus_y_minus) >= 0;
C2 = abs(Dy_minus_y_minus) <= abs(Dy_plus_y_minus);
C = Dy_minus + delta_y_half * (C1.*C2 .* Dy_minus_y_minus+ (C1&~C2) .* Dy_plus_y_minus);

D1 = (Dy_plus_y_plus .*  Dy_plus_y_minus) >= 0;
D2 = abs(Dy_plus_y_plus) <= abs(Dy_plus_y_minus);
D = Dy_plus  - delta_y_half * (D1.*D2 .* Dy_plus_y_plus  + (D1&~D2) .* Dy_plus_y_minus);
        
Amax = A > 0;
Amin = A < 0; 
Bmax = B > 0;
Bmin = B < 0;
Cmax = C > 0;
Cmin = C < 0;
Dmax = D > 0;
Dmin = D < 0;

grad_x = Amax .* A + Bmin .* B;  
grad_y = Cmax .* A + Dmin .* B; 

delta_plus  =  sqrt((Amax.*A).^2 + (Bmin.*B).^2 + (Cmax.*C).^2 + (Dmin.*D).^2);
delta_minus =  sqrt((Bmax.*B).^2 + (Amin.*A).^2 + (Dmax.*D).^2 + (Cmin.*C).^2);


% for i=1:i_end
%     for j=1:j_end
%         if (Dx_minus_x_minus(i,j)  * Dx_plus_x_minus(i,j)) < 0
%             m = 0;
%         elseif abs(Dx_minus_x_minus(i,j)) <= abs(Dx_plus_x_minus(i,j))
%             m = Dx_minus_x_minus(i,j);
%         else
%             m = Dx_plus_x_minus(i,j);
%         end
%         A = Dx_minus(i,j) + delta_x_half * m;
%         
%         
%         if (Dx_plus_x_plus(i,j) * Dx_plus_x_minus(i,j)) < 0
%             m = 0;
%         elseif abs(Dx_plus_x_plus(i,j)) <= abs(Dx_plus_x_minus(i,j))
%             m = Dx_plus_x_plus(i,j);
%         else
%             m = Dx_plus_x_minus(i,j);
%         end
%         B = Dx_plus(i,j)  - delta_x_half * m;
%         
%         if (Dy_minus_y_minus(i,j) * Dy_plus_y_minus(i,j)) < 0
%             m = 0;
%         elseif abs(Dy_minus_y_minus(i,j)) <= abs(Dy_plus_y_minus(i,j))
%             m = Dy_minus_y_minus(i,j);
%         else
%             m = Dy_plus_y_minus(i,j);
%         end
%         C = Dy_minus(i,j) + delta_y_half * m;
%         
%         if (Dy_plus_y_plus(i,j) * Dy_plus_y_minus(i,j)) < 0
%             m = 0;
%         elseif abs(Dy_plus_y_plus(i,j)) <= abs(Dy_plus_y_minus(i,j))
%             m = Dy_plus_y_plus(i,j);
%         else
%             m = Dy_plus_y_minus(i,j);
%         end
%         D = Dy_plus(i,j)  - delta_y_half * m;
%        
%         
%         
%         grad_x(i,j) = max(A,0) + min(B,0);
%         grad_y(i,j) = max(C,0) + min(D,0);
% 
%         delta_plus(i,j)  =  sqrt(max(A,0)^2 + min(B,0)^2+...
%             max(C,0)^2 + min(D,0)^2);
% 
%         delta_minus(i,j) =  sqrt(max(B,0)^2 + min(A,0)^2+...
%             max(D,0)^2 + min(C,0)^2);
% 
% 
% 
% %       grad_x(i,j) = max(A(i,j),0) + min(B(i,j),0);
% %       grad_y(i,j) = max(C(i,j),0) + min(D(i,j),0);
% %       
% %       delta_plus(i,j)  =  sqrt(max(A(i,j),0)^2 + min(B(i,j),0)^2+...
% %                                max(C(i,j),0)^2 + min(D(i,j),0)^2);
% %       
% %       delta_minus(i,j) =  sqrt(max(B(i,j),0)^2 + min(A(i,j),0)^2+...
% %                                max(D(i,j),0)^2 + min(C(i,j),0)^2);    
% 
% 
%    end
% end

% diff = sum(sum(delta_plus_r - delta_plus));
% diff2 = sum(sum(delta_minus_r - delta_minus));
% p=1;
% w1=[1 -1];
% w2=[1 -2 1];
% 
% Dy_minus          = filter(w1,delta_y,phi,[],1); % you have to define first element
% Dy_minus(1,:)     = Dy_minus(2,:);
% 
% Dy_plus           = Dy_minus; % add last element!!! 
% Dy_plus(1,:)      = [];
% Dy_plus(end+1,:)  = Dy_plus(end,:);
% 
% Dy_minus_y_minus  =  filter(w2,delta_y*delta_y,phi,[],1); % you have to define first two element
% Dy_minus_y_minus(1,:) = Dy_minus_y_minus(3,:);
% Dy_minus_y_minus(2,:) = Dy_minus_y_minus(3,:);
% 
% Dy_plus_y_minus   =  Dy_minus_y_minus; % add last element!!! 
% Dy_plus_y_minus(1,:)=[];
% Dy_plus_y_minus(end+1,:) = Dy_plus_y_minus(end,:);
% 
% Dy_plus_y_plus    =  Dy_plus_y_minus; %  add last element!!! 
% Dy_plus_y_plus(1,:)=[]; 
% Dy_plus_y_plus(end+1,:) = Dy_plus_y_plus(end,:);
% 
% 
% 
% Dx_minus          =  filter(w1,delta_x,phi,[],2); % xou have to define first element
% Dx_minus(:,1) = Dx_minus(:,2);
% 
% Dx_plus           =  Dx_minus; % add last element!!! 
% Dx_plus(:,1)=[];
% Dx_plus(:,end+1)=Dx_plus(:,end);
% Dx_minus_x_minus  =  filter(w2,delta_x*delta_x,phi,[],2); % xou have to define first two element
% Dx_minus_x_minus(:,1) = Dx_minus_x_minus(:,3);
% Dx_minus_x_minus(:,2) = Dx_minus_x_minus(:,3);
% 
% Dx_plus_x_minus   =  Dx_minus_x_minus; % add last element!!! 
% Dx_plus_x_minus(:,1)=[]; 
% Dx_plus_x_minus(:,end+1) = Dx_plus_x_minus(:,end);
% 
% Dx_plus_x_plus    =  Dx_plus_x_minus; %  add last element!!! 
% Dx_plus_x_plus(:,1)=[];
% Dx_plus_x_plus(:,end+1) = Dx_plus_x_plus(:,end);
%
%
%
% function m = switch_m(D1, D2)
% 
% if (D1 * D2) < 0
%    m = 0; 
% elseif abs(D1) <= abs(D2)
%    m = D1;
% else
%    m = D2;
% end