function[corfac]=circumferenceCorrFactor(xx,yy,rr,msx,msy)
%circumference correction calculates a vector containing the correction factor
%(for edge correction in Ripley's k-function) for values of rr
%circumference correction: fraction of circumference of circle centered at
%point P=(xx,yy) with radius rr (inside rectangular image) falling into the 
%rectangle - this fraction becomes smaller as the point gets closer to one
%of the rectangle's edges, and as the radius of the circle increases
%if the circle falls completely inside the rectangle, the value is zero

%Notes: - this function assumes that rr is a vector
%       - this version of the function allows rr to be larger than ms/2

% SYNOPSIS   [corfac]=circumferenceCorrFactor(xx,yy,rr,msx,msy)
%
% INPUT:    xx      x-position 
%           yy      y-position
%           rr      distance (vector)
%           msx     image size in x
%           msy     image size in y
%           
% OUPUT:    corfac  correction factor; vector of same length as rr %       
%
% Dinah Loerke, January 27, 2007


% while xx and yy are the point coordinates, x and y are the *distances* from
% the nearest edge of the image in that direction; e.g. x=xx for a point close
% to the left edge, and x=(msx-xx) for a point close to the right edge 
x=min(xx,(msx-xx));
y=min(yy,(msy-yy));

%xo and yo are, conversely, the maximum distance to the image edge
xo=max(xx,(msx-xx));
yo=max(yy,(msy-yy));

rmax=max(size(rr));
%correction factor is initialized
corfac=ones(rmax,1);

for i=1:rmax
    r=rr(i);
    
    %the distance r needs to be positive, and smaller than the maximum
    %distance of the current point (xx,yy) from any of the edges, since
    %else no fraction of circumference falls inside the image. for larger 
    %distances, the correction factor is technically zero, which creates a
    %continuous function for increasing values of r, but this function
    %should not generate zero-value outputs, since the main function 
    %divides a defined value by the correction factor, so that we'd get 
    %a 'divide by zero' error message. thus, the default value for this 
    %situation is set to nan instead of zero
    if ( (r>0) & (r<sqrt(xo^2+yo^2)) )
             
            %consider contributions of all 4 quadrants separately
            %each quadrant can contribute at most 0.25, at least 0 to the
            %total corrfactor
            cont=zeros(4,1);
            % loop over all four quadrants
            for qq=1:4
                switch qq
                    %determine relevant distances from the edges (xd,yd) 
                    %for each quadrant
                    case(1)
                        % quadrant 1 : upper left 
                        xd=xx;
                        yd=msy-yy;
                    case(2)
                        % quadrant 2 : upper right 
                        xd=msx-xx;
                        yd=msy-yy;
                    case(3)
                        % quadrant 3 : lower right 
                        xd=msx-xx;
                        yd=yy;
                    case(4)
                         % quadrant 4 : lower left 
                        xd=xx;
                        yd=yy;
                end % of switch
                
                %now xd and yd are the distances from the closest image 
                %edges in the current quadrant
                %if both xd and yd are smaller than r, this means that the 
                %circle goes across both edges of the image
                if((xd<r)&&(yd<r))
            
                    %if additionally r is smaller than the distance to the corner,
                    %we have case 1: the circle cuts off the corner
                    if( sqrt(xd^2+yd^2) > r )
                        %two angles on the side of the corner that define
                        %the fraction of the circumference outside the
                        %image
                        out_angle1 = acos(xd/r);
                        out_angle2 = acos(yd/r);
                        %inside angle is 90-outside angles
                        in_angle = (pi/2)-(out_angle1+out_angle2);
                        %contribution is 0.25*relative fraction of inside 
                        %angle in this quadrant
                        cont(qq)=0.25*(in_angle/(pi/2));
                    
                    %else the circle clears the corner (case 2); in this
                    %case, no part of the circumference of this quadrant
                    %lies inside the image, and the contribution is zero
                    else
                        cont(qq)=0;
                    end % of if
                    
                %else, either x or y is larger than r (or both are)
                %this means that the circle goes across at most one edge of
                %the image in this quadrant (case 3)
                else
                    % z is whatever is the smaller distance to the edge; if
                    % the quarter circle is fully inside the image in this 
                    % quadrant, then z=r                
                    z= min(min(xd,yd),r) ;
                    cont(qq)=0.25 * ( asin(z/r) / (pi/2) );
                end % of if-else
            end % of for qq
            
            %add contributions of four quadrants
            corfac(i)=sum(cont);
    %else corfac is zero (e.g. for point identity, because center point is
    %not counted), or undefined because the distance r is not within the 
    %limits of the image for this particular point (xx,yy) 
    else
        if (r==0)
            corfac(i)=1;
        else
            corfac(i)=nan;
        end
    end % of if
end % of for i


end  % of function
