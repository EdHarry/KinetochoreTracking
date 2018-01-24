function[cpar,pvr,dpvr, cpar2]=clusterQuantRipley(mpm,imsizex,imsizey,norm,ra);
% clusterQuantRipley calculates a quantitative clustering parameter based on
% Ripley's K-function (a spatial statistics function for point patterns)
%
% SYNOPSIS   [cpar,pvr,dpvr]=clusterQuantRipley(mpm,imsizex,imsizey);
%       
% INPUT      mpm:   mpm file containing (x,y) coordinates of points in the
%                   image in succesive columns for different time points
%            imsizex:   x-size of the image (maximum possible value for x-coordinate)
%            imsizey:   y-size of the image (maximum possible value for
%                       y-coordinate)
%            norm:   number of images for normalization; e.g. number of
%                   images before adding growth factor to cells
%            ra: number of planes for rolling average - the function will
%                return cpar vectors that have the same size (i.e. the same
%                number of time points) as the original mpm, the missing
%                edge points will be filled up with duplicates
%
%            NOTE: in Johan's mpm-files, the image size is 1344 x 1024
%               pixels
%            NOTE2: this function uses Ripley's circumference correction
%            NOTE3: Although this function is not actually named after 
%                   Lt. Ellen Ripley, she certainly would deserve to have 
%                   a kick-ass matlab function named after her.
%
%
% OUTPUT     cpar:  for each time point, a single clustering parameter 
%                   value is extracted from the pvr function  
%            pvr:   for each plane (time point), the function calculates
%                   the function pvr=points vs radius, i.e. the number of
%                   points contained in a circle of increasing radius
%                   around an object, averaged over all objects in the
%                   image
%                   NOTE: the default size for the radius implicit in the 
%                   pvr function is [1,2,3...,minimsize] where minimsize is the
%                   smalller dimension of imsizex,imsizey
%                   This function corresponds to Ripley's K-function (/pi)
%            dpvr:  H(r) function
%            cpar2: additional parameters
%
% DEPENDENCES   ClusterQuantRipley uses {pointsincircle,clusterpara}
%               ClusterQuantRipley is used by { }
%
% Dinah Loerke, October 7th, 2004
% last update: April 05,2005

%if not specified, ra = 1, so that no rolling average is performed
%for values >1, the function always uses the closest odd number for rolling
%average, i.e. 3,5,7,...
if nargin < 5, ra = 1; end
ra=1+2*ceil((ra-1)/2);

%create vector containing x- and y-image size
matsiz=[imsizex imsizey];
%rs is size of distance vector in Ripley function
rs=round(min(matsiz)/2);

%determine size of mpm-file
[nx,ny]=size(mpm);
%number of frames
numframes = round(ny/2);

%initialize results matrix pvr; x-dimension equals the employed number of
%values for the circle radius, y-dimension equals number of planes of the
%input mpm-file
pvr=zeros(rs,numframes);
dpvr=zeros(rs,numframes);
cpar=[1:numframes];
cpar2=[1:numframes];

%initialize temporary coordinate matrix matt, which contains the object 
%coordinates for one plane of the mpm
matt=zeros(nx,2);


%%======================================================================
%
%    CALCULATE (relative) NUMBER OF POINTS IN CIRCLE OF INCREASING SiZE
%
%=======================================================================


%cycle over all planes of series, using two consecutive columns of mpm input
%matrix as (x,y) coordinates of all measured points
for k=1:(round(ny/2))
    
    %matt is set to two consecutive columns of input matrix m1
    matt(:,:)=mpm(:,(2*k-1):(2*k));
    
    %since the original mpm file contains a lot of zeros, these zeros are 
    %deleted in the temporary coordinate matrix to yield a matrix containing
    %only the nonzero points of matt, smatt
    [nz1,e]=size(nonzeros(matt(:,1)));
    [nz2,e]=size(nonzeros(matt(:,2)));
    % dovar is do-variable to determine whether function is performed on
    % this plane or not (due to missing objects or non-matching coordinates)
    dovar=1;
    if( (nz1==nz2) && (nz1>0) )
        smatt=[nonzeros(matt(:,1)), nonzeros(matt(:,2)) ];
    else
        dovar=0;
        disp(['Error in frame ',num2str(k), ' of input mpm']);
        if(nz1~=nz2)
            disp(['unequal number of entries for x and y-coordinates']);
        end
        if (nz1==0) 
        disp(['no objects (i.e. no nonzero entries) in this frame']);
        end
    end
    
    %comment/uncomment the next five lines if you want to monitor progress
    %prints number of objects for every 10th line
    [smx,smy]=size(smatt);
    tempnp=max([smx,smy]);
    if(mod(k,10)==0)
        disp(['plane ',num2str(k),'   number of objects ', num2str(tempnp)]);
    end  % of if

    if (dovar>0)
    %now determine number of objects in circle of increasing radius,
    %averaged over all objects in smatt, and normalized with point density
    %tempnp/(msx*msy)
        [pvrt,nump]=pointsincircle(smatt,matsiz);
    %result is already normalized with point density tempnp/(msx*msy)
        pvr(:,k)=pvrt(:);
    
    end  % of if
        
end %of for

%%======================================================================
%
%    CALCULATE H(r) FUNCTION AND CLUSTER PARAMETERS 
%
%=======================================================================

%From the calculated function pvrt (number of points versus circle
%radius), we calculate a quantitative clustering parameters, cpar.
%If specified by user input as ra>1, use rolling average of pvr function
%In this case, we fill the missing start and end points with duplicates -
%since some frames get lost by the rolling average, but we want the
%output to have the same number of frames as the original, so that in
%subsequent analysis, we don't have to correct the normpoint (e.g. the time
%point of drug addition) for the rolling average magnitude.
%The "padding" in front and back (variable "shift") is 1 each for ra=3, 
%2 each for ra=5, etc.
%For example, for ra=5, the averages for the first 3 frames are all 1:5, then 2:6,
% etc. For the last frames, the averages are n-5:n-1, then n-4:n in the last
% three frames.
shift=round((ra-1)/2);
maxp=round(ny/2);
for k=1:numframes
    %startpoint cannot be below zero, is minumum of k-shift
    astartpoint=max(1,k-shift);
    %endpoint is shifted to right from startpoint, but cannot exceed 
    %the maximum length of the vector
    aendpoint=min(maxp,astartpoint+ra-1);
    %startpoint is restrosepctively adjusted to endpoint
    astartpoint=aendpoint-(ra-1);
    points=[k astartpoint aendpoint];
    %disp(points);
    %calculate rolling average of pvrt
    pvrt = mean( pvr(:,astartpoint:aendpoint),2 );
    smatt=[nonzeros(matt(:,1)), nonzeros(matt(:,2)) ];
    [smx,smy]=size(smatt);
    tempnp=max([smx,smy]);
    [cpar(k),dpvr(:,k), cpar2(k)]=clusterpara(pvrt,k);
    
end % of for

%normalize cpar with initial value, either first point or specified lebgth
%of normalization norm before drug addition
if(norm>length(cpar))
    norm=1;
end
normfac=nanmean(cpar(1:norm));
cpar=cpar/normfac;

fpnormval=nanmean(cpar2(1:norm));
cpar2=cpar2/fpnormval;
%for the second cluster parameter, firstpoint, exclude points that are
%significantly larger than twice the starting value
cpar2(cpar2>2.4)=nan;


end


function[cpar1, dpvrt, cpar2]=clusterpara(pvrt,k);
%clusterpara calculates a quantitative cluster parameter from the input
%function (points in circle) vs (circle radius)
% SYNOPSIS   [cpar]=clusterpara(pvrt);
%       
% INPUT      pvrt:   function containing normalized point density in 
%                   circle around object
%                   spacing of points implicitly assumes radii of 1,2,3...
%            k= number of plane in series (is used for monitoring progress)
%
% OUTPUT     cpar:    cluster parameter
%            dpvrt:   difference function of p vs r
%
% DEPENDENCES   clusterpara uses {DiffFuncParas}
%               clusterpara is used by {FractClusterQuant}
%
% Dinah Loerke, September 13th, 2004

%============================================================
% calculate difference L(d)-d function, using L(d)=sqrt(K(d)),
% since K(d) is already divided by pi
% This difference function is from now on called H(r)
%=============================================================
len=max(size(pvrt));
de=(1:len);
Hr=pvrt-de'.^2;

% for difffuncparas, we want to extract inclination around central point;
% for high degree of cell division, the normalization of Hr affects this
% inclination; therefore, to conserve the height of the first rise 
% corresponding to the close neighborhood, we scale with number of cells; 
% 
%numc=max(nump);
%dpvrt=Hr*tempnp^(1.1);
dpvrt=Hr;

%extract parameters from diff
[cpar1,cpar2]=DiffFuncParas(dpvrt,k);

end    

    
function[p1,p2]=DiffFuncParas(Hr,k)
%DiffFuncParas calculates a number of quantitative cluster parameter 
%from the input function, the difference function
% SYNOPSIS   DiffFuncParas(diff);
%       
% INPUT      Hr:  difference function as calculated in clusterpara
%                   vector with len number of points
%%          k:  plane in series
% OUTPUT     p1,p2:  cluster parameters
%                    currently: p1=inclination of the first rise of the
%                                   H(t) function (clustering) 
%                               of total clustering)
%                               p2= position of first rise
%
% DEPENDENCES   DiffFuncParas uses {}
%               DiffFuncParas is used by {clusterpara}
%
% Dinah Loerke, September 13th, 2004


%% firstpoint: point where diff function systematically rises above zero 
%% definition: diff>0 AND (diff)'>0 to exlude noisy one-point rises above
%% zero. if there exists no such point (for completely scattered
%% distributions), firstpoint is set to nan and incl is set to zero


vec=Hr;

% smooth with Gaussian filter
xs=-8:1:8;
amps=exp(-(xs.^2)/(2*(2.5^2)));
namps=amps/sum(amps);
filtershape = namps;
%the above definition of the filter will introduce a shift to the filtered
%vector, this will need to be compensated by a counter-shift firther below
shift=8;
[filtervec] = filter(filtershape,1,vec);


% the filter shifts the function several points to the right, this is
% compensated by removing first points (number of points is=shift)
% the size of the filter excludes variation caused by small cell numbers,
% e.g. the wedge artifact for single cell increases at small distances, if
% the number of lonely frames isn't too large
filtervec(1:shift)=[];

%devc=first differential
dvec=diff(filtervec);
%ddvec = second differential
ddvec=diff(dvec);

% detvec= determination vector; has the function to differentiate between 
% clustered distributions, where we can calculate cluster parameters, and 
% scattered distrubutions where this is impossible. detvec is set to zero 
% where either the original function is below zero (indicating scattering)
% or where the inclination (of the filtered function) is below zero (as 
% would be the case for a single isolated above-zero data point, which is
% not followed shortly after by an additional data point - the 'shortly
% after' depends on the range of the filtering

%inclination of curve has to be > 0...
detvec1=dvec;
detvec1(1:20)=0;
detvec1(dvec<0)=0;
minimum=min(find(detvec1));

%...and the value of the function has to be > 0, too.
detvec3=filtervec;
detvec3(1:minimum)=0;
detvec3(filtervec<0)=0;
firstzerocrosspoint=min(find(detvec3));

%for determination of inclination:
detvec2=ddvec;
detvec2(1:minimum)=0;
detvec2(ddvec>0)=0;

if (length(nonzeros(detvec3))>1)
    firstpoint=firstzerocrosspoint;
    % beginning point begp and end point endp for calculating inclination of the
    % function diff
    begp=firstpoint-2;
    endp=begp+round(firstpoint/2);
    if ( endp>length(dvec) )
        endp=length(dvec);
    end
    inclination=mean(dvec(begp:endp));
    
%% =====================================================================
%% if anything is funny with the results of the analysis,
%  uncomment the following paragrpah for a display of the single traces 
%  during determination
%% ===========================================================
     plot(vec,'b.');
     axis([ 20 120 -5000 15000]);
     hold on
     plot(filtervec,'b-');
     ypts=[vec(begp) vec(endp)];
     xpts=[begp endp];
     plot(xpts, ypts, 'r.');
     plot(firstpoint,filtervec(firstpoint),'go');
     text(30,5000,num2str(k));
     pause(0.1);
     hold off;
    
else
    firstpoint=NaN;
    inclination=0;
end

%disp(['firstpoint ',num2str(firstpoint)]);
p1=inclination;
p2=firstpoint;

end



function[m2,num]=pointsincircle(m1,ms)
%pointsincircle calculates the average number of points in a circle around
%a given point as a function of the circle radius (averaged over all points
%and normalized by total point density); this function is called Ripley's
%K-function in statistics, and is an indication of the amount of clustering
%in the point distribution
% 
% SYNOPSIS   [m2]=pointsincircle(m1,ms);
%       
% INPUT      m1:   matrix of size (n x 2) containing the (x,y)-coordinates of n
%                  points
%            ms: vector containing the parameters [imsizex imsizey] (the 
%                   x-size and y-size of the image)
%            NOTE: in Johan's mpm-files, the image size is 1344 x 1024
%               pixels
%
%
% OUTPUT     m2:    vector containing the number of points in a circle 
%                   around each point, for an increasing radius;
%                   radius default values are 1,2,3,....,min(ms)
%                   function is averaged over all objects in the image
%           num: number of points
%
% DEPENDENCES   pointsincircle uses {distanceMatrix, circumferenceCorrectionFactor}
%                   (distanceMatrix, circumferenceCorrectionFactor added to this file)
%               pointsincircle is used by {FractClusterQuant }
%
% Dinah Loerke, October 4th, 2004


[lm,wm]=size(m1);

%for points at the edges (where the circle of increasing size is cut off by
%the edges of the image), this function corrects for the reduced size of 
%the circle using the 
%function circumferenceCorrectionFactor

msx=ms(1);
msy=ms(2);
minms=min(ms);
rs=round(minms/2);

%create neighbour matrix m3
%matrix m3 contains the distance of all points in m1 from all points
%in itself
[mdist]=distanceMatrix(m1,m1);

%create numpoints vector (number of points in circle of corresponding radius)
%loop over all radius values between 1 and minms
%initialize m2 vector
%m2=1:rs;

%create corrections factor matrix (same dimension as mdist)
%contains correction factor for precise radii (point distances) around each
%point; for the zero entry at identity (p11,p22,p33), cfm equals one
corrFacMat=ones(lm);
for n=1:lm
    corrFacMat(n,:) = circumferenceCorrectionFactor(m1(n,1),m1(n,2),mdist(n,:),msx,msy);
end
thresh_mdist = mdist;

for r=rs:-1:1
    %for given radius, set all values of mdist higher than the radius value
    %to zero
    
    thresh_mdist(thresh_mdist > r) = 0;
    
    %count all leftover points equally => set to one
    
    thresh_mdistones = thresh_mdist;
    thresh_mdistones(thresh_mdist > 0) = 1;
    %what's num represent, again...
    num(r)=sum(thresh_mdistones(:))/2;  
    %weight every counted point with the circumference correction factor
    %calculated previously in corrFacMat
    tempfinal=thresh_mdistones./corrFacMat;    
    %sum over entire matrix to get number of points
    npv=sum(tempfinal(:));
      
    %to average, divide sum by number of points (=columns)
    npv=(npv/lm);
    
    %in order to be able to quantitatively compare the clustering in 
    %distributions of different point densities, this npv value must now 
    %be corrected for overall point density, which is lm/msx*msy; the 
    %resulting normalized function is (if we also divide by pi to scale for
    %the circle area) more or less a simple square function;
    %it is a perfect square function for a perfectly random distribution of
    %points
    m2(r)=npv/(pi*(lm-1)/(msx*msy));
    %using (lm-1) and not lm is Marcon&Puech's correction (2003)
    
end

end
  

function[m2]=distanceMatrix(c1,c2)
%this subfunction makes a neighbour-distance matrix for input matrix m1
%input: c1 (n1 x 2 points) and c2 (n2 x 2 points) matrices containing 
%the x,y coordinates of n1 or n2 points
%output: m2 (n1 x n2) matrix containing the distances of each point in c1 
%from each point in c2
[ncx1,ncy1]=size(c1);
[ncx2,ncy2]=size(c2);
m2=zeros(ncx1,ncx2);
for k=1:ncx1
    for n=1:ncx2
        d=sqrt((c1(k,1)-c2(n,1))^2+(c1(k,2)-c2(n,2))^2);
        m2(k,n)=d;
    end
end
end
    
function[corfac]=circumferenceCorrectionFactor(xx,yy,rr,msx,msy)
%circumference correction calculates a vector containing the correction factor
%(for edge correction in Ripley's k-function) for values of rr
%circumference correction: fraction of circumference of circle centered at
%point P=(xx,yy) with radius rr (inside rectangular image) falling into the 
%rectangle - this fraction becomes smaller as the point gets closer to one
%of the rectangle's edges, and as the radius of the circle increases
%if the circle falls completely inside the rectangle, the value is zero

%1. this function assumes that rr is a vector

% SYNOPSIS   [corfac]=circumferenceCorrectionFactor2(xx,yy,rr,msx,msy)
%       
% Dinah Loerke, October 6, 2004


x=min(xx,(msx-xx));
y=min(yy,(msy-yy));

rmax=max(size(rr));
corfac=ones(rmax,1);

for i=1:rmax
    r=rr(i);
    %if both x and y are smaller than r
    if(r>0)
        
        if((x<r)&&(y<r))
            if( (x<r) && (y<r) && (sqrt(x^2+y^2)>r) )
                corfac(i)=(2*asin(x/r)+2*asin(y/r))/(2*pi);
            else
                corfac(i)=(0.5*pi+asin(x/r)+asin(y/r))/(2*pi);
            end
         %if either x OR y OR neither is smaller than r 
        else
             z=min( min(x,y),r );
            corfac(i)=(pi+2*asin(z./r))/(2*pi);
        end
    end
end

end

