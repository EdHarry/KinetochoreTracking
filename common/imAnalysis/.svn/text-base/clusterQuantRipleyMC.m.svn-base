function[cpar1,cpar2, cpar3,pvr,dpvr]=clusterQuantRipleyMC(mpm,imsizex,imsizey,norm,ra,corr);
% clusterQuantRipleyMC calculates a quantitative clustering parameter based on
% Ripley's K-function (a spatial statistics function for point patterns)
% This function corrects for population growth (by cell division) by using
% a Monte-Carlo simulation of the culster growth
% SYNOPSIS   [cpar1,cpar2,
% cpar3,pvr,dpvr]=clusterQuantRipleyMC(mpm,imsizex,imsizey,norm,ra,corr);
%       
% INPUT      mpm:       mpm file containing (x,y) coordinates of points in the
%                       image in succesive columns for different time points
%            imsizex:   x-size of the image (maximum possible value for x-coordinate)
%            imsizey:   y-size of the image (maximum possible value for
%                       y-coordinate)
%            norm:      number of images for normalization; e.g. number of
%                       images before adding growth factor to cells
%            ra:        number of planes for rolling average - the function will
%                       return cpar vectors that have the same size (i.e. the same
%                       number of time points) as the original mpm, the missing
%                       edge points will be filled up with ; a typical
%                       value is e.g. ra=5
%           corr:       optional input, enter 'flip' for overedge flip 
%                       correction, the default is ripley's correction
%                       (using correction factor determined by the fraction
%                       of the circumference in the image)
%
%            NOTE:      in Johan's mpm-files, the image size is 1344 x 1024
%                       pixels
%            NOTE2:     this function uses ripley's correction unless
%                       otherwise specified
%            NOTE3:     Although this function is not actually named after 
%                       Lt. Ellen Ripley, she certainly would deserve to have 
%                       a kick-ass matlab function named after her.
%
%
% OUTPUT    
% for each time point, a single clustering parameter value is extracted
% from the pvr function
%                   
%           cpar1:  total integrated value of dpvr (H(r)) from 0 to rs
%                   (value of rs is specified below, depends on image size) 
%           cpar2:  first point for dpvr to cross from negative to positive
%                   values - means the distance for which a significantly
%                   increased number of neighbors can be found in
%                   comparison to statistical distributions
%           cpar3:  partial integrated value of dpvr (H(r)) from 0 to rs
%           pvr:    for each plane (time point), the function calculates
%                   the function pvr=points vs radius, i.e. the number of
%                   points contained in a circle of increasing radius
%                   around an object, averaged over all objects in the
%                   image
%                   
%                   NOTE: the default size for the maximum radius in the
%                   function has been changed!! The maximum size used to be
%                   minsiz/2, i.e. half of the smaller dimension of the
%                   image; in the current version with an updated function
%                   for edge correction, the default value for max r is
%                   half of the diagonal length.
%                   This function corresponds to Ripley's K-function (/pi)
%
%           dpvr:   H(r) function, defined as Hr=pvrt-de'.^2 (in the
%                   literature, you will also find a variant of this
%                   function usually called the L(d)-d function, which is 
%                   L(d)=sqrt(pvrt)-de'
%            
%
% DEPENDENCIES      clusterQuantRipleyMC uses the external functions:
%                   - nearNeiDistProbHisto
%                   - numPointsinMpm
%                   - mpmMCaddsim
%                                                
%
% Dinah Loerke, last update: Feb 07, 2005





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SET DEFAULT VALUES AND INITIALIZE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if not otherwise specified, ra = 1, so that no rolling average is performed
% for values >1, the function always uses the closest odd number for rolling
% average, i.e. 3,5,7,...
if nargin < 5, ra = 1; end
ra=1+2*ceil((ra-1)/2);

if ra>norm
    disp('potential source of error: moving average window (ra) is larger than normalization window (norm)!');
end

% create vector containing x- and y-image size
matsiz=[imsizex imsizey];

%number of cycles for MC simulation
simLoop = 3;

% rs is size of distance vector in Ripley function
% in former versions, rs was half the smaller image side length, in this
% new version, rs is half the length of the diagonal - which means that for
% all distances up to rs, at least some fraction of the circle lies within
% the image, i.e. it is possible to detect neighbors at distance r even for
% the centermost points in the image.
% if this version is changed to include even larger distances, the
% pointsincircle function also has to be changed, because then the case may
% occur that for points in the center of the image, the entire circle lies
% outside the image and no new points may be detected

%rs=round(sqrt((imsizex*imsizey)/pi));
global rs
rs = round(sqrt(imsizex^2+imsizey^2)/2);

% determine size of mpm-file
[nx,ny]=size(mpm);
if ((ny/2) < norm)
    error('ERROR - normalization window is larger than number of frames in input mpm!');
end

% number of frames 
numframes = round(ny/2);
numcells = numPointsinMpm(mpm);

% initialize results matrix pvr; x-dimension equals the employed number of
% values for the circle radius, y-dimension equals number of planes of the
% input mpm-file
pvr=zeros(rs,numframes);
pvr_growthCorr=zeros(rs,numframes);

dpvr=zeros(rs,numframes);
dpvr_growthCorr=zeros(rs,numframes);

cpar=[1:numframes];
cpar2=[1:numframes];
cpar3=[1:numframes];

cpar_growthCorr=[1:numframes];
cpar2_growthCorr=[1:numframes];
cpar3_growthCorr=[1:numframes];

% initialize temporary coordinate matrix matt, which contains the object 
%coordinates for one plane of the mpm
matt=zeros(nx,2);

corrVar = 0;
if nargin > 5
    if (corr == 'flip')
        corrVar = 1;
        disp('edge correction: toroidal (flip image over edge)');
    else
        disp('edge correction: Ripley (weigh with ratio of circumference inside image)');
    end
end


%%======================================================================
%  
%   CALCULATE K(r)= Ripley's K-function of measured distribution
%   (number of points in circle of radius r)
%
%=======================================================================

disp('calculating Ripley K-function for input mpm...');
[pvr]=RipleyKfunc(mpm,matsiz,corrVar);



%======================================================================
%
%    For growth correction, simulate a growing (but non-scattering and
%    non-moving) cell population 
%    construct MPM and then also calculate Ripley's K-function   
%
%=======================================================================

np = numcells;
numaddcells=zeros(length(np)-1,1);

% modify np (remove point disappearances if necessary) and make addcell 
% vector for MC simulation
for k=1:length(np)
    if (k>1)
        if ( np(k)<np(k-1)) 
            np(k)=np(k-1);
        end
        numaddcells(k-1) = np(k)-np(k-1);
    end
end

% determine probability distribution for Monte Carlo simulation from as many
% frames as possible, i.e. for the number of frames specified
% for normalization

pdavmat=zeros(100,norm);
for k=1:norm
    tx=2*k-1;
    ty=2*k;
    mpmktemp=mpm(:,tx:ty);
    [pd,pdfilter]=nearNeiDistProbHisto(mpmktemp);
    pdavmat(:,k)=pdfilter;
end
pdav=mean(pdavmat,2);
% the approximate cell diameter is the distance for which the probability
% distribution has its maximum value
cellDiam = min(find(pdav ==(max(pdav))));


% make simulated mpm using average probablity distribution pdav
% use original data for frames 1-norm
% start simulating @ frame 1
% start simulating @ frame randomly chosen from numbers 1:norm-1

% average over 3 simulated mpms
pvr_growthCorrMat=zeros(rs,numframes,3);

for simnum=1:simLoop
    
    startnum = norm-1;
    %startnum = ceil((norm-1)*rand(1));
    mpmstart=mpm(:,(2*startnum-1):(2*startnum));
    mpmcurr=mpmstart;
    mpmtot=zeros(nx,2*numframes);
    mpmtot(:,1:2*norm)=mpm(:,1:2*norm);

    for k=(1+startnum):numframes
        [mpmMC]=mpmMCaddsim(mpmcurr,[imsizex imsizey], numaddcells(k-1),pdav);
        [xmc,ymc]=size(mpmMC);
        mpmtot(1:xmc,(2*k-1):(2*k))=mpmMC; 
        mpmcurr = mpmMC;
    end % of for

    %now calculate Ripley's K-function for simulated mpm just as for the original

    disp(['calculating Ripley K-function for simulated mpm number ',num2str(simnum)]);
    [pvr_growthCorrTemp]=RipleyKfunc(mpmtot,matsiz,corrVar);
    pvr_growthCorrMat(:,:,simnum) = pvr_growthCorrTemp;
    
end % of for simnum
pvr_growthCorr=mean(pvr_growthCorrMat,3);
% pvr_growthCorr is now the averaged pvr function for a non-scattering and
% non-moving dividing cell population



%%======================================================================
%
%    CALCULATE H(r) FUNCTION AND CLUSTER PARAMETERS 
%
%=======================================================================

%From the calculated function pvr (number of points versus circle
%radius), we calculate a number of quantitative clustering parameters, cpar.
%
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
% interpoint distance for homogeneous distribution of this density
%ipdist=round(sqrt(((imsizex*imsizey)./numcells)/0.86));


%now do final analysis including filtering (which uses forward and backward
%data, so it has do be done in a separate loop
h = waitbar(0,'H(r) calculation');

figure
for k=1:numframes
    waitbar(k/numframes);
    %disp(['final parameter determination in frame ',num2str(k)]);
    %startpoint cannot be below zero, is minumum of k-shift
    astartpoint=max(1,k-shift);
    %endpoint is shifted to right from startpoint, but cannot exceed 
    %the maximum length of the vector
    aendpoint=min(maxp,astartpoint+ra-1);
    %startpoint is restrosepctively adjusted to endpoint
    astartpoint=aendpoint-(ra-1);
    points=[k astartpoint aendpoint];

    %calculate rolling average of pvrt
    pvrt = mean( pvr(:,astartpoint:aendpoint),2 );
    pvrt_growthCorr = mean( pvr_growthCorr(:,astartpoint:aendpoint),2 );
    
    len=max(size(pvrt));
    de=(1:len);
%============================================================
% calculate difference function H(r)=K(r)-r^2
% difference L(d)-d function, using L(d)=sqrt(K(d)),
% since K(d) is already divided by pi
% This difference function is from now on called H(r)
%In the context of this function, it's called dpvr
%=============================================================
    
    dpvr(:,k) = pvrt - de'.^2;
    switch k
        case 1
            plot(dpvr(:,k),'b-');
            hold on
        case norm
            plot(dpvr(:,k),'c-');
        case norm+round((numframes-norm)/2)
            plot(dpvr(:,k),'g-');
        case numframes
            plot(dpvr(:,k),'r-');
    end % of switch
            
    dpvr_growthCorr(:,k) = pvrt_growthCorr - de'.^2;
    
end
close(h);

%use average of 1:norm to determine border values for paremetr
%determination if desired

dpvr_rest = mean(dpvr(:,1:norm),2);
drest_smooth = dpvr_rest;
flen = round(cellDiam/4);
for f = flen+1:(len-flen); drest_smooth(f) = mean(dpvr_rest(f-flen:f+flen)); end

d_drest = diff(drest_smooth);
dd_drest = diff(d_drest);

d2(1:2*cellDiam) = 0;
rangef = drest_smooth(2*cellDiam:10*cellDiam);
if ( min(rangef)<0 )
    restCSiz = 2*cellDiam + min(find(rangef<0));
else
    restCSiz = 2*cellDiam + find(rangef == min(rangef));
end
%restCSiz is the first point for which the function drops below zero again,
%i.e. something like the maximum cluster size at rest

%hwb = waitbar(0,'parameter calculation');

    
for k=1:numframes
    %waitbar(k/numframes);
    tempnp=numcells(k);
    % calculate cluster parameters in subfunctions for original and
    % simulated pvr
    %subplot(1,2,1);
    [cpar1(k),cpar2(k),cpar3(k)]=clusterpara(dpvr(:,k),matsiz,cellDiam,restCSiz);
    %subplot(1,2,2);
    [cpar1_growthCorr(k),cpar2_growthCorr(k),cpar3_growthCorr(k)]=clusterpara(dpvr_growthCorr(:,k),matsiz,cellDiam,restCSiz);
    
end % of for
%close(hwb);
%h2 = figure;

%=========================================================================
% 
%  Normalize and correct measured values using the normalization phase and
%  using the results of the growth correction analysis
%
%==========================================================================

%normalize cpar with initial value, either first point or specified lebgth
%of normalization norm before drug addition
if(norm>length(cpar))
    norm=1;
end

%c1 is total integrated area

normfac1=nanmean(cpar1(1:norm));
cpar1=cpar1/normfac1;
normfac1_growthCorr=nanmean(cpar1_growthCorr(1:norm));
cpar1_growthCorr=cpar1_growthCorr/normfac1_growthCorr;

%c2 ist first point of rise above zero

normfac2=nanmean(cpar2(1:norm));
cpar2=cpar2/normfac2;
normfac2_growthCorr=nanmean(cpar2_growthCorr(1:norm));
cpar2_growthCorr=cpar2_growthCorr/normfac2_growthCorr;


%cpar 3 ist partial integrated area

normfac3=nanmean(cpar3(1:norm));
cpar3=cpar3/normfac3;
normfac3_growthCorr=nanmean(cpar3_growthCorr(1:norm));
cpar3_growthCorr=cpar3_growthCorr/normfac3_growthCorr;


%===============================
%
%       Growth correction
%
%===============================

% to correct integrals for the value of growth correction, smooth 
% smooth growthCorr functions
% since the simulated vectors correspond to a non-scattering distribution, 
% they don't go below zero and just reflect the decay of the initial height 
% of the H(r) function due to normalization effects
% Thus, the measured value is divided by the (normalized) simulation value


%filter (smooth) and shift to account for filter effects 
filterLength = ceil(norm/4);
preVec3 = ones(1,filterLength);
preVec3(:) = cpar3_growthCorr(1);
cpar3_growthCorrSmooth = filter( ones(1,filterLength)/filterLength,1,[preVec3 cpar3_growthCorr] );
cpar3_growthCorrSmooth(1:round(1.5*filterLength)-1)=[];
l3 = length(cpar3_growthCorrSmooth);
cpar3_growthCorrSmooth(l3+1:numframes)=cpar3_growthCorrSmooth(l3);

%same for preVec1 - this procedure is not necessary for cpar2, since this
%parameter is intrinsically independent of cluster growth
preVec1 = ones(1,filterLength);
preVec1(:) = cpar1_growthCorr(1);
[cpar1_growthCorrSmooth] = filter( ones(1,filterLength)/filterLength,1,[preVec1 cpar1_growthCorr] );
cpar1_growthCorrSmooth(1:round(1.5*filterLength)-1)=[];
l1 = length(cpar1_growthCorrSmooth);
cpar1_growthCorrSmooth(l1+1:numframes)=cpar1_growthCorrSmooth(l1);

%uncomment the following paragraphs to display results, don't forget the
%trailing line 421

figure
plot(cpar3,'ro');
axis([0 numframes -0.5 1.5]);
hold on
plot(cpar3_growthCorr,'b-');

cpar3=cpar3./cpar3_growthCorrSmooth;
cpar1=cpar1./cpar1_growthCorrSmooth;

plot(cpar3,'g.');



%step 2: for the second cluster parameter, firstpoint, exclude those points
%that occur for strongly scattered distributions, i.e. where the value of
%cpar3 have dropped below a set threshold of e.g. 10%, since in the very
%scattered distributions, this parameter eventually becomes meaningless and
%very noisy
%
cpar2(find(cpar3<0.1))=nan;


end % main function





%==========================================================================

function[cpar1,cpar2,cpar3]=clusterpara(Hr,matsiz,cellDiam,restCSiz);
%clusterpara calculates a quantitative cluster parameter from the input
%function (points in circle) vs (circle radius)
% SYNOPSIS   [cpar1,cpar2,cpar3]=clusterpara(Hr,matsiz,cellDiam,restCSiz);
%       
% INPUT      Hr:    function containing H(r) normalized point density in 
%                   circle around object for specific t
%                   spacing of points implicitly assumes radii of 1,2,3...
%            matsiz: size of image
%            cellDiam = cell diamater (from nearest neighbor distance
%            probability histogram)
%            restCSiz : resting whole cluster size   
%               
%
% OUTPUT     cpar1:  total integral
%            cpar2:  position of first rise from <0 to >0
%            cpar3:  partial integral from ~cellDiam to size of first 
%            cluster (which is defined from at rest)
%
% DEPENDENCIES   clusterpara uses {}
%               
%
% Dinah Loerke, August 20th 2005


%% firstpoint: point where diff function systematically rises above zero 
%% definition: diff>0 AND (diff)'>0 to exlude noisy one-point rises above
%% zero. if there exists no such point (for completely scattered
%% distributions), firstpoint is set to nan and incl is set to zero

global rs

vec = Hr;
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

%since in this version of the function, the growth is already corrected
%for, the ipdist value here can remain fixed
%for ipdist, we choose a value related to the cell diameter, so that cells
%of different size can be compared to each other
ipdist=min (5*cellDiam, restCSiz);
restZeroCross = round(0.85 * cellDiam);
intlength = round(0.6*(ipdist-restZeroCross));
entIntegral = sum(abs(Hr(restZeroCross:rs)));
parIntegral = sum(Hr(restZeroCross:restZeroCross+intlength));    


% if (length(nonzeros(detvec3))>1)
%     firstpoint=firstzerocrosspoint;
%     % beginning point begp and end point endp for calculating inclination of the
%     % function diff
%     begp=firstpoint-2;
%     endp=begp+round(firstpoint/2);
%     if ( endp>length(dvec) )
%         endp=length(dvec);
%     end
%     inclination=mean(dvec(begp:endp));
%     %parIntegral = sum(Hr(firstpoint:ipdist));
%     
% %% =====================================================================
% %% if anything is funny with the results of the analysis,
% %  uncomment the following paragrpah for a display of the single traces 
% %  during determination
% %% ===========================================================
%      
%           
%      plot(filtervec,'b-');
%      axis([ 0 ipdist+20 -2000 17000]);
%      hold on
%      ypts=[filtervec(cellDiam) filtervec(ipdist)];
%      xpts=[cellDiam ipdist];
%      plot(xpts, ypts, 'r.');
%      ypts2=[filtervec(restZeroCross)];
%      xpts2=[restZeroCross];
%      plot(xpts2, ypts2, 'g.');
%      
%      pause(0.05);
%      hold off;
%     
% else
%     firstpoint=NaN;
%     inclination=0;
% end


if (parIntegral>0)
    firstpoint = firstzerocrosspoint;
%    parIntegral = sum(Hr(firstzerocrosspoint:firstzerocrosspoint+intlength));    
else
    firstpoint = NaN;
end


%cpar=inclination;
cpar1 = entIntegral;
cpar2 = firstpoint;
cpar3 = parIntegral;


end





%==========================================================================


function[m2]=pointsincircle(m1,ms,corrVar)
%pointsincircle calculates the average number of points in a circle around
%a given point as a function of the circle radius (averaged over all points
%and normalized by total point density); this function is called Ripley's
%K-function in statistics, and is an indication of the amount of clustering
%in the point distribution
% 
% SYNOPSIS   [m2]=pointsincircle(m1,ms,corrVar)
%       
% INPUT      m1:   matrix of size (n x 2) containing the (x,y)-coordinates of n
%                  points
%            ms: vector containing the parameters [imsizex imsizey] (the 
%                   x-size and y-size of the image)
%            corr: 'flip' or 'rip' correction factor
%
%            NOTE: in Johan's mpm-files, the image size is 1344 x 1024
%               pixels
%
%
% OUTPUT     m2:    vector containing the number of points in a circle 
%                   around each point, for an increasing radius;
%                   radius default values are 1,2,3,....,min(ms)
%                   function is averaged over all objects in the image
%           
%
% DEPENDENCES   pointsincircle uses {distanceMatrix, circumferenceCorrectionFactor}
%                   (distanceMatrix, circumferenceCorrectionFactor added to this file)
%               pointsincircle is used by {FractClusterQuant }
%
% Dinah Loerke, October 4th, 2004


[lm,wm]=size(m1);

%for points at the edges (where the circle of increasing size is cut off by
%the edges of the image), this function corrects for the reduced size of 
%the circle by assuming a contuous distribution, i.e. by flipping the image
%over the edge

msx=ms(1);
msy=ms(2);

global rs

%=========== for toroidal correction

if (corrVar==1)
    
    midxo = round(0.6*msx);
    midxu = round(0.4*msx);
    midyo = round(0.6*msy);
    midyu = round(0.4*msy);

    %duplicate of m1 with all surrounding fields included
    mq5=m1;
    pos = find ((m1(:,1)>midxu) & (m1(:,2)>midyu));
    mq1=[(m1(pos,1)-msx) (m1(pos,2)-msy)];
    pos = find (m1(:,2)>midyu);
    mq2=[(m1(pos,1))     (m1(pos,2)-msy)];
    pos = find ((m1(:,1)<midxo) & (m1(:,2)>midyu));
    mq3=[(m1(pos,1)+msx) (m1(pos,2)-msy)];
    pos = find (m1(:,1)>midxu);
    mq4=[(m1(pos,1)-msx) (m1(pos,2))     ];
    pos = find (m1(:,1)<midxo);
    mq6=[(m1(pos,1)+msx) (m1(pos,2))     ];
    pos = find ((m1(:,1)>midxu) & (m1(:,2)<midyo));
    mq7=[(m1(pos,1)-msx) (m1(pos,2)+msy)];
    pos = find (m1(:,2)<midyo);
    mq8=[(m1(:,1))     (m1(:,2)+msy)];
    pos = find ((m1(:,1)<midxo) & (m1(:,2)<midyo));
    mq9=[(m1(:,1)+msx) (m1(:,2)+msy)];

    mq=[mq1; mq2; mq3; mq4; mq5; mq6; mq7; mq8; mq9];
    %disp( ['matrix = ',num2str(length(m1)),'  mq = ',num2str(length(mq))]);
    %create neighbour matrix m3
    %matrix mdist contains the distance of all points in m1 from all points
    %in itself
    [mdist]=distanceMatrix(m1,mq);
    

    % ===================== else ripley correction

else  %default is 'rip'=correction 'cause it's faster
    %create neighbour matrix m3
    %matrix m3 contains the distance of all points in m1 from all points
    %in itself
    [mdist]=distanceMatrix(m1,m1);

    % %create corrections factor matrix (same dimension as mdist)
    % %contains correction factor for precise radii (point distances) around each
    % %point; for the zero entry at identity (p11,p22,p33), cfm equals one
    corrFacMat=ones(lm);
    for n=1:lm
        corrFacMat(n,:) = circumferenceCorrectionFactor(m1(n,1),m1(n,2),mdist(n,:),msx,msy);
    end
    
    
end %of if-else

thresh_mdist = mdist;

for r=rs:-1:1
    %for given radius, set all values of mdist higher than the radius value
    %to zero
    
    thresh_mdist(thresh_mdist > r) = 0;
    
    %count all leftover points equally => set to one
    
    thresh_mdistones = thresh_mdist;
    thresh_mdistones(thresh_mdist > 0) = 1;
    
    %weight every counted point with the circumference correction factor
    %calculated previously in corrFacMat
    if (corrVar == 0)
        tempfinal=thresh_mdistones./corrFacMat;
    else
        tempfinal=thresh_mdistones;  
    end
           
    %sum over entire matrix to get number of points
    npv = nansum(tempfinal(:));
    
    
    %to average, divide sum by number of points (=columns)
    %identity points have value 1 in the correction factor matrix, but they
    %are NOT counted in the thresholded matrix
    npv=npv/lm;
    
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

end % of function
  




%==========================================================================

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




%==========================================================================

function[corfac]=circumferenceCorrectionFactor(xx,yy,rr,msx,msy)
%circumference correction calculates a vector containing the correction factor
%(for edge correction in Ripley's k-function) for values of rr
%circumference correction: fraction of circumference of circle centered at
%point P=(xx,yy) with radius rr (inside rectangular image) falling into the 
%rectangle - this fraction becomes smaller as the point gets closer to one
%of the rectangle's edges, and as the radius of the circle increases
%if the circle falls completely inside the rectangle, the value is zero

%Notes: - this function assumes that rr is a vector
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





%==========================================================================

function[pvr]=RipleyKfunc(mpm,imsiz,corrFac);
%[pvr]=RipleyKfuncFlip(mpm,imsiz); 
%pure Kfunction calculation from mpm



%create vector containing x- and y-image size
imsizex = imsiz(1);
imsizey = imsiz(2);

%rs is size of distance vector in Ripley function, defined in main function
global rs

%determine size of mpm-file
[nx,ny]=size(mpm);
%number of frames
numframes = round(ny/2);
%numcells = zeros(numframes,1);

%initialize results matrix pvr; x-dimension equals the employed number of
%values for the circle radius, y-dimension equals number of planes of the
%input mpm-file
pvr=zeros(rs,numframes);

%initialize temporary coordinate matrix matt, which contains the object 
%coordinates for one plane of the mpm
matt=zeros(nx,2);

%%======================================================================
%
%    CALCULATE (relative) NUMBER OF POINTS IN CIRCLE OF INCREASING SiZE
%
%=======================================================================

h = waitbar(0,'Ripley (K-function) calculation');
%cycle over all planes of series, using two consecutive columns of mpm input
%matrix as (x,y) coordinates of all measured points
for k=1:(round(ny/2))
    waitbar(k/round(ny/2));
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
        disp(['Error in frame ',num2str(k), ' of input mpm:']);
        if(nz1~=nz2)
            error(['unequal number of entries for x and y-coordinates']);
        end
        if (nz1==0) 
        error(['no objects (i.e. no nonzero entries) in this frame']);
        end
    end
    
    %comment/uncomment the next five lines if you want to monitor progress
    %prints number of objects for every 10th line
%     [smx,smy]=size(smatt);
%     numcells(k)=max([smx,smy]);
%     if(mod(k,10)==0)
%         disp(['plane ',num2str(k),'   number of objects ', num2str(numcells(k))]);
%     end  % of if

    if (dovar>0)
    %now determine number of objects in circle of increasing radius,
    %averaged over all objects in smatt, and normalized with point density
    %tempnp/(msx*msy)
        [pvrt]=pointsincircle(smatt,imsiz,corrFac);
    %result is already normalized with point density tempnp/(msx*msy)
        pvr(:,k)=pvrt(:);
    
    end  % of if
        
end %of for
close(h);

end % of function




