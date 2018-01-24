function [kr,lr]=RipleysKfunction(mpm1,mpm2,imsiz,dist,corrFacMat,normArea);
% RipleysKfunction calculates Ripley's K-function for a given MPM,
% allowing cross-corrlation between two MPMs
% SYNOPSIS  [kr,lr]=RipleysKfunction(mpm,imsiz,dist,corrFacMat);
%       
% INPUT     mpm1:      mpm file containing (x,y) coordinates of points in 
%                      the image in succesive columns for different time 
%                      points
%           mpm2:      second mpm
%           imsiz:     x,y-size of the image (maximum possible value for x-coordinate)
%           dist:      distance vector e.g. [1:50]
%           corrMat:   OPTIONAL if a correction matrix is pre-calculated 
%                      outside of this function, it can be used directly
%           normArea:  OPTIONAL if the point density for normalization
%                      should not be based on number of points per total
%                      rectangular image size, but e.g. per a different
%                      area of interest (a mask of which may have already
%                      been used to calculate corrMat), then the area
%                      for normalization should be entered here
%                                 
%           NOTE: IF you want to pre-calculate the correction factor matrix
%           (which is recommended because it saves time), then run the 
%           function: 
%           [corrFacMatrix] = makeCorrFactorMatrix(imsiz, dist, samplesize,
%           mask); 
%           
%           IF nargin == 5              the provided corrFacMat is used
%           IF corrFacMat is an integer corrFacMat is calculated using
%                                       this integer value as samplesize
%           IF nargin == 4              Ripley's edge correction is used
%           IF nargin == 3 (or ==nan)   dist is set as default [1:rs]           
%
% OUTPUT                    
%           kr:     Ripley's K-function 
%           lr:     Besag's L-function = sqrt(K(r))-r
%           for both, every columns contains the kr/lr function for one
%           frame of the mpm-file
%                                                
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: dloerke
% 
% last modified
% DATE: 29-Jan-2008 (last update)
%
%

% create vector containing x- and y-image size
imsizex = imsiz(1);
imsizey = imsiz(2);

% if no distance vector is chosen explicitly or if dist is nan, the 
% distance vector is set as default to 1:rs, where rs is the half diagonal 
% of the image (this is the standard in the literature)
rs = round(sqrt(imsizex^2+imsizey^2));
if nargin<4
    distvec = [1:rs];
    nr = rs;
elseif isnan(dist)
    distvec = [1:rs];
    nr = rs;
else
    distvec = dist;
    nr = length(distvec);
end

% if corrFacMat is specified but consists only of an integer value, 
% calculate the matrix using that value as sample, else set corrFacMat to
% empty
if (nargin<5)
    corrFacMat = [];
else
    if length(corrFacMat)==1
        sample = corrFacMat;
        corrFacMat = makeCorrFactorMatrix(imsiz,sample);
    end
end

%determine size of mpm-file
[nx1,ny1]=size(mpm1);
%number of frames
numframes1 = round(ny1/2);
%determine size of mpm-file
[nx2,ny2]=size(mpm2);
%number of frames
numframes2 = round(ny2/2);

% function requires both mpms to have the same number of frames
if numframes1~=numframes2
    error('number of frames in 2 mpms doesn''t match');
else
    numf = numframes1;
end

%initialize results matrix pvr; x-dimension equals the employed number of
%values for the circle radius, y-dimension equals number of planes of the
%input mpm-file
kr = zeros(nr,numf);
lr = kr;

% loop over all frames
for i=1:numf
    % relevant nonzero x,y-positions in mpm1
    pos1 = find(mpm1(:,i*2)>0);
    cmpm1 = mpm1(pos1,2*i-1:2*i);
    np1 = length(pos1);
    
    % relevant nonzero x,y-positions in mpm2
    pos2 = find(mpm2(:,i*2)>0);
    cmpm2 = mpm2(pos2,2*i-1:2*i);
    np2 = length(pos2);
    
    fprintf(' frame %04d',i);
    
    % if there are any relevant points in this frame - at least one point
    % each for cross-correlation (different matrices mpm1/mpm2), at least 
    % two points for auto-corr (mpm1 == mpm2)
    if ( min(np1,np2)>0 ) & ( (np1+np2)>2 )
        
        % output pr: # of points as a function of distance    
        [pr,nump] = pointsincircleCross(cmpm1,cmpm2,imsiz,distvec,corrFacMat);
    
        % normalized pr - normalize by total point density
        totaldensity = (nump)/(imsizex*imsizey); 
        if nargin>5
            totaldensity = nump/normArea;
        end
    
        prnorm = (pr/pi)*(1/totaldensity);
        kr(:,i) = prnorm;
        lr(:,i) = sqrt(prnorm) - distvec;
    else
        kr(:,i) = nan*dist;
        lr(:,i) = nan*dist;
    end
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b');
    
end % of for i-loop

fprintf('\n');

end % of function





%=========================================================================
%=========================================================================
%=========================================================================
%====================       SUBFUNCTIONS    ==============================
%=========================================================================
%=========================================================================
%=========================================================================




function [npvr,nump]=pointsincircleCross(m1,m2,ms,dist,corrMat)
% pointsincircle calculates the average number of points in a circle around
% a given point as a function of the circle radius (averaged over all points
% and normalized by total point density); this function is called Ripley's
% K-function in statistics, and is an indication of the amount of clustering
% in the point distribution
% 
% SYNOPSIS   [m2,num]=pointsincircle(m1,m2,ms,dist,corrMat);
%       
% INPUT     m1:     matrix of size (n1 x 2) containing the (x,y)-coordinates 
%                   of n1 points; these points are considered the CHILDREN
%           m2:     matrix of size (n2 x 2) containing the (x,y)-coordinates 
%                   of n2 points; these points are considered the PARENTS
%           ms:     vector containing the parameters [imsizex imsizey] (the 
%                   x-size and y-size of the image)
%           dist:   distance vector
%           corrMat:   correctionFactor matrix; if this matrix is empty or
%                   nargin<5, then the simple Ripley's edge correction is
%                   used to calculate the correction factor
%
%
% OUTPUT    npvr:   vector containing the number of points in a circle 
%                   around each point, for an increasing radius;
%                   radius default values are 1,2,3,....,min(ms)
%                   function is averaged over all objects in the image
%           nump:   number of points
%
% Dinah Loerke, Jan 29th, 2008

[nump1,numd1]=size(m1);
[nump2,numd2]=size(m2);

msx=ms(1);
msy=ms(2);

% in the following, the matrix mdist will contain the distance of all 
% points in m1 from all points in m2, where m1 are the children and m2 are
% the parents; because of the way the function below is set up in terms of 
% rows-columns, the order in the DistanceMatrix function needs to be 
% parent-child

[mdist]=DistanceMatrix(m2,m1);

% NOTE: In the old version of the Ripley, since it was designed for
% self-correlation, only distances > 0 were considered in the distance
% histogram below. In this version, since it allows cross-correlation,
% zero-distances have to be included for DIFFERENT mpms, but should be
% excluded for IDENTICAL mpms

% identity variable 
ivar = (m2==m1);

% monitor progress
%fprintf(' progress %02d',0);

if isempty(m1) | isempty(m2) | isempty(mdist) | (mdist==0)
    npvr = nan*dist;
    nump = 0;
else
    
    for i=1:length(mdist(:,1))

        % matrix histmat contains the distance histogram for each point, where 
        % every row represents the distance histogram for the cell at that 
        % position
        histmat(i,:) = histc(nonzeros(mdist(i,:)),[0 dist]);

        % determine the circumference correction vector for this point; in the 
        % matrix corrmat, the columns represent the points, and the rows 
        % represent the entries for the radius vector dist
        % IF corrMat exists, use these values directly, else calculate corrmat 
        % in situ with Ripley's edge correction
        if ~isempty(corrMat)
            cpx = round(m2(i,1)); cpy = round(m2(i,2));
            corrmat(i,:)= corrMat(cpx,cpy,:);
        else
            % x- and y-distances of this center point (the parent point) from the
            % nearest edge 
            ex = min(m2(i,1),1+msx-m2(i,1));
            ey = min(m2(i,2),1+msy-m2(i,2));
            corrmat(i,:)= circumferenceCorrectionFactor(ex,ey,dist,msx,msy);
        end

        % update iter 
        iter = round( 100*(i/length(mdist(:,1))) );
        %if iter<100, fprintf('\b\b%02d',iter); end

    end



    % multiply the matrices histmat and the correction vector
    histmat(:,length(dist)+1)=[];
    pointsmat = histmat./corrmat;

    % for every point, the (number-of-points within radius r)-function is the
    % cumulative sum of the values in the corresponding row
    pointsincirclemat = cumsum(pointsmat,2);

    % average for this frame over all existing points
    npvr = nanmean(pointsincirclemat,1);
    nump = nump1;

end % of if there are any points



end % of function
  


%==========================================================================

function [corfac]=circumferenceCorrectionFactor(xx,yy,rr,msx,msy)
% circumference correction calculates a vector containing the correction factor
% (for edge correction in Ripley's k-function) for values of rr
% circumference correction: fraction of circumference of circle centered at
% point P=(xx,yy) with radius rr (inside rectangular image) falling into the 
% rectangle - this fraction becomes smaller as the point gets closer to one
% of the rectangle's edges, and as the radius of the circle increases
% if the circle falls completely inside the rectangle, the value is zero

% Notes: - this function assumes that rr is a vector
%       - this version of the function allows rr to be larger than ms/2

% SYNOPSIS   [corfac]=circumferenceCorrectionFactor(xx,yy,rr,msx,msy)
%       
% Dinah Loerke, January 27, 2006


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


%%======================================================================





function [m2]=DistanceMatrix(c1,c2);
%this subfunction makes a neighbour-distance matrix for input matrix c1
%(n1 x 2 points) and c2
%output: m2 (n1 x n1) matrix containing the distances of each point in c1 
%from each point in c2

[np1,sd1]=size(c1);
[np2,sd2]=size(c2);

m2=zeros(np1,np2);

for k = 1:np1
    for n = 1:np2
        d = sqrt((c1(k,1)-c2(n,1))^2+(c1(k,2)-c2(n,2))^2);
        m2(k,n)=d;
    end
end


end

