function [trackClass,mssSlope,genDiffCoef,scalingPower,normDiffCoef] ...
    = trackMSSAnalysis(tracks,probDim,momentOrders,alphaMSS)
%TRACKMSSANALYSIS classifies trajectories based on their moment scaling spectrum

%SYNPOSIS [trackClass,mssSlope,genDiffCoef,scalingPower,normDiffCoef] ...
%    = trackMSSAnalysis(tracks,probDim,momentOrders,alphaMSS)
%
%INPUT  tracks      : Matrix indicating the positions and amplitudes of the
%                     tracked features. Number of rows = number of tracks,
%                     number of columns = 8*number of time points.
%                     Each row consists of
%                     [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...].
%                     NaN is used to indicate time points where the track
%                     does not exist.
%       probDim     : Problem dimensionality. Optional. Default: 2.
%       momentOrders: Orders of moments to be calculated.
%                     Optional. Default: 0 through 6.
%       alphaMSS    : Alpha-value for classification. Can take the values
%                     0.2, 0.1, 0.05 and 0.01.
%                     Optional. Default: 0.1.
%
%OUTPUT trackClass  : # tracks x 1 vector of track classification.
%                     Values mean the following ...
%                     0 = stalled. (NOT IMPLEMENTED YET)
%                     1 = confined Brownian.
%                     2 = pure Brownian.
%                     3 = Brownian with drift (directed).
%                     NaN = not classified.
%       mssSlope    : # tracks x 1 vector of each track's slope of the line
%                     representing moment scaling power vs. moment order.
%                     NaN indicates tracks that could not be analyzed.
%       genDiffCoef : # tracks x # orders array of generalized diffusion
%                     coefficients for every moment order considered.
%                     NaN indicates tracks that could not be analyzed.
%       scalingPower: # tracks x # orders array of powers with which moment
%                     values scale with time.
%                     NaN indicates tracks that could not be analyzed.
%       normDiffCoef: # tracks x 1 vector of each track's "normal"
%                     diffusion coefficient.
%                     NaN indicates tracks that could not be analyzed.
%
%REMARKS (1) Algorithm is based on Ewers et al. 2005. PNAS 102:
%15110-15115 and Ferrari et al. 2001. Physica D 154: 111-137.
%(2) Analysis assumes that there are no kinks in the moment scaling
%spectrum curve, i.e. that the motion is strongly self-similar. Weakly
%self-similar processes will generate an MSS which is piece-wise
%continuous, hence before fitting to estimate the slope the curve must be
%chopped into smaller straight-line pieces (but this is not done).
%
%Khuloud Jaqaman, March 2008

%% input

if nargin < 1
    disp('--trackMSSAnalysis: Please input tracks to analyze.');
    return
end

if nargin < 2 || isempty(probDim)
    probDim = 2;
end

if nargin < 3 || isempty(momentOrders)
    momentOrders = 0 : 6;
end
numOrders = length(momentOrders);

if nargin < 4 || isempty(alphaMSS)
    alphaMSS = 0.1;
end

%get number of tracks and frames
[numTracks,numFramesMovie] = size(tracks);
numFramesMovie = numFramesMovie / 8;

%find indices of tracks that are >= 20 frames long - do not attempt
%to calculate moments for shorter tracks
criteria.lifeTime.min = 20;
indx4diff = chooseTracks(tracks,criteria);
clear criteria

%% alpha-value for classification

%determine threshold based on alpha-value and dimensionality
switch probDim
    case 1
        switch alphaMSS
            case 0.2 %10th percentile and 90th percentile
                [mssThreshNeg,mssThreshPos] = threshMSS1D_p20(numFramesMovie);
            case 0.1 %5th percentile and 95th percentile
                [mssThreshNeg,mssThreshPos] = threshMSS1D_p10(numFramesMovie);
            case 0.05 %2.5th percentile and 97.5th percentile
                [mssThreshNeg,mssThreshPos] = threshMSS1D_p05(numFramesMovie);
            case 0.01 %0.5th percentile and 99.5th percentile
                [mssThreshNeg,mssThreshPos] = threshMSS1D_p01(numFramesMovie);
        end
    case 2
        switch alphaMSS
            case 0.2 %10th percentile and 90th percentile
                [mssThreshNeg,mssThreshPos] = threshMSS2D_p20(numFramesMovie);
            case 0.1 %5th percentile and 95th percentile
                [mssThreshNeg,mssThreshPos] = threshMSS2D_p10(numFramesMovie);
            case 0.05 %2.5th percentile and 97.5th percentile
                [mssThreshNeg,mssThreshPos] = threshMSS2D_p05(numFramesMovie);
            case 0.01 %0.5th percentile and 99.5th percentile
                [mssThreshNeg,mssThreshPos] = threshMSS2D_p01(numFramesMovie);
        end
    case 3
        switch alphaMSS
            case 0.2 %10th percentile and 90th percentile
                [mssThreshNeg,mssThreshPos] = threshMSS3D_p20(numFramesMovie);
            case 0.1 %5th percentile and 95th percentile
                [mssThreshNeg,mssThreshPos] = threshMSS3D_p10(numFramesMovie);
            case 0.05 %2.5th percentile and 97.5th percentile
                [mssThreshNeg,mssThreshPos] = threshMSS3D_p05(numFramesMovie);
            case 0.01 %0.5th percentile and 99.5th percentile
                [mssThreshNeg,mssThreshPos] = threshMSS3D_p01(numFramesMovie);
        end
end

%% memory for trajectory classification

%classification means ...
%0 = stalled
%1 = confined Brownian
%2 = pure Brownian
%3 = drift/directed
%NaN = unclassified

trackClass = NaN(numTracks,1);
mssSlope = NaN(numTracks,1);
genDiffCoef = NaN(numTracks,numOrders);
scalingPower = NaN(numTracks,numOrders);
normDiffCoef = NaN(numTracks,1);

%% moments and their scaling with time

for iTrack = indx4diff'

    %get track start and end time
    trackSEL = getTrackSEL(tracks(iTrack,:));
    startTime = trackSEL(1);
    endTime = trackSEL(2);
    numTimePoints = trackSEL(3);

    %extract track's coordinates and their standard deviations
    coordinates = [tracks(iTrack,1:8:end)' tracks(iTrack,2:8:end)' tracks(iTrack,3:8:end)'];
    coordinates = coordinates(startTime:endTime,:);
    standardDevs = [tracks(iTrack,5:8:end)' tracks(iTrack,6:8:end)' tracks(iTrack,7:8:end)'];
    standardDevs = standardDevs(startTime:endTime,:);

    %define maximum time lag for moment calculation
    maxLag = min(30,floor(numTimePoints/4));

    %calculate track moments
    trackMomentsT = calcTrackMoments(coordinates,standardDevs,momentOrders,maxLag);
    trackMoments = [trackMomentsT.momentValues];
    trackMoments = trackMoments(:,1:2:end);

    %estimate the moment scaling spectrum (MSS),
    %i.e. the scaling power for all moments
    scalingPowerT = NaN(1,numOrders);
    genDiffCoefT = NaN(1,numOrders);
    for iOrder = 1 : length(momentOrders)

        %caculate ln(lag) and ln(moment)
        lnTime = log((1:maxLag)');
        lnMoment = log(trackMoments(:,iOrder));
        
        %remove any NaNs
        indxGood = find(~isnan(lnMoment));
        lnTime = lnTime(indxGood);
        lnMoment = lnMoment(indxGood);
        
        %if there are moments to fit ...
        if length(lnMoment) > 1

            %fit a straight line in the plot of lnMoment vs. lnTime
            slParam = polyfit(lnTime,lnMoment,1);

            %get scaling power and generalized diffusion coefficient
            scalingPowerT(iOrder) = slParam(1);
            genDiffCoefT(iOrder) = exp(slParam(2)) / 2 / probDim;
            
            %if this is the 2nd moment, calculate the "normal" diffusion
            %coefficient
            if momentOrders(iOrder)==2
                options = optimset('Display','off','Jacobian','on');
                lnSlope = lsqcurvefit(@strLineFun2,1,lnTime(1:min(5,...
                    length(lnTime))),lnMoment(1:min(5,length(lnMoment))),...
                    [],[],options);
                normDiffCoefT = exp(lnSlope) / 2 / probDim;                
            end
            
        end

    end

    %keep only non-NaN scaling powers
    indxGood = find(~isnan(scalingPowerT));
    momentOrders4fit = momentOrders(indxGood);
    scalingPowerT = scalingPowerT(indxGood);
    genDiffCoefT = genDiffCoefT(indxGood);
    
    %if there are non-NaN scaling powers
    if ~isempty(scalingPowerT)

        %fit a straight line to the MSS
        slParam = polyfit(momentOrders4fit,scalingPowerT,1);

        %get the slope of the line
        mssSlopeT = slParam(1);

        %classify track as ...
        %1 = confined Brownian, if MSS slope < mssThreshNeg
        %2 = pure Brownian, if mssThreshNeg <= MSS slope <= mssThreshPos
        %3 = directed, if MSS slope > mssThreshPos
        if ~isnan(mssSlopeT)
            if mssSlopeT < mssThreshNeg(numTimePoints)
                trackClass(iTrack) = 1;
            elseif mssSlopeT > mssThreshPos(numTimePoints)
                trackClass(iTrack) = 3;
            else
                trackClass(iTrack) = 2;
            end
        end

        %save additional output information
        mssSlope(iTrack) = mssSlopeT;
        genDiffCoef(iTrack,:) = genDiffCoefT;
        scalingPower(iTrack,:) = scalingPowerT;
        normDiffCoef(iTrack) = normDiffCoefT;
        
    end

end

%% subfunction 1
function [y,d] = strLineFun2(logSlope,x)

y = logSlope + x;
d = ones(size(x));


%% thresholds

function [mssThreshNeg,mssThreshPos] = threshMSS1D_p20(nTP)

%1D, alpha = 0.2

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 60 100 500];
slopeM = [0.002434794446981 0.000542795021531 0.000165329365168 0];
slopeP = [-0.001541239220102 -0.00036860800289 -0.00002638219148 0];
interseptM = [0.130233710469096 0.243753675996077 0.319246807268674 0.40191148985285];
interseptP = [0.672499246818346 0.602141373785626 0.567918792644692 0.554727696904506];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS1D_p10(nTP)

%1D, alpha = 0.1

%threshold curve parameters
turnPointsM = [20 50 200 500];
turnPointsP = [20 60 150 500];
slopeM = [0.00403865723972 0.00073054296121 0.000167789786604 0];
slopeP = [-0.001863373027824 -0.000207607837099 -0.000068295975124 0];
interseptM = [0.018637170190375 0.184042884115918 0.296593519037011 0.380488412339042];
interseptP = [0.727469044719516 0.628123133276013 0.607226353979748 0.5730783664177];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS1D_p05(nTP)

%1D, alpha = 0.05

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 50 200 500];
slopeM = [0.004221335524384 0.000794743062117 0.000192611669839 0];
slopeP = [-0.002478698219027 -0.000203422663306 -0.000101280672122 0];
interseptM = [-0.059561585354062 0.146033962381945 0.266460240837543 0.362766075757204];
interseptP = [0.772420366381382 0.658656588595335 0.638228190358568 0.587587854297634];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS1D_p01(nTP)

%1D, alpha = 0.01

%threshold curve parameters
turnPointsM = [20 50 150 500];
turnPointsP = [20 45 100 500];
slopeM = [0.00748519566052008 0.00132812830129705 0.000291477165467116 0];
slopeP = [-0.00263804313196125 -0.000652822904393821 -0.000123920022338378 0];
interseptM = [-0.27508360480994 0.0327697631512118 0.188267433525701 0.334006016259259];
interseptP = [0.839817228313045 0.750482318072511 0.697592029866967 0.635632018697778];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS2D_p20(nTP)

%2D, alpha = 0.2

%threshold curve parameters
turnPointsM = [20 50 200 500];
turnPointsP = [20 60 150 500];
slopeM = [0.002365010936411 0.000327173962804 0.000137961052668 0];
slopeP = [-0.001158501020215 -0.000170757301414 -0.000033415212947 0];
interseptM = [0.227688919307209 0.329580767987585 0.367423350014689 0.436403876348735];
interseptP = [0.634413262884249 0.575148639756167 0.554547326486174 0.537839720012491];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS2D_p10(nTP)

%2D, alpha = 0.1

%threshold curve parameters
turnPointsM = [20 50 200 500];
turnPointsP = [20 50 100 500];
slopeM = [0.002834926966999 0.00040315841383 0.000171843106892 0];
slopeP = [-0.001737592454918 -0.000337743432604 -0.000073809631405 0];
interseptM = [0.169732430047339 0.291320857705762 0.337583919093415 0.423505472539345];
interseptP = [0.684919949284652 0.614927498168961 0.588534118049036 0.551629302346427];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS2D_p05(nTP)

%2D, alpha = 0.05

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 50 100 500];
slopeM = [0.002405756738527 0.000449948051639 0.000201793444222 0];
slopeP = [-0.001937063703976 -0.00055853262239 -0.000095409685751 0];
interseptM = [0.14034289465302 0.257691415866306 0.307322337349635 0.40821905946088];
interseptP = [0.724427540398996 0.6555009863197 0.609188692655777 0.561483849780384];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS2D_p01(nTP)

%2D, alpha = 0.01

%threshold curve parameters
turnPointsM = [20 50 150 500];
turnPointsP = [20 60 200 500];
slopeM = [0.00435441385708182 0.000812259003939792 0.000273699080289468 0];
slopeP = [-0.00172943051183377 -0.000284948443095389 -0.000153857046768591 0];
interseptM = [-0.00562462883719271 0.171483113819909 0.252267102367457 0.389116642512191];
interseptP = [0.779260212544801 0.692591288420498 0.666373009155138 0.589444485770843];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS3D_p20(nTP)

%3D, alpha = 0.2

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 60 100 500];
slopeM = [0.001359327197525 0.000280732353905 0.000108640351462 0];
slopeP = [-0.000746906632969 -0.000348357342053 -0.000037074457479 0];
interseptM = [0.299405762325546 0.364121452942741 0.39853985343145 0.452860029162242];
interseptP = [0.609183824504497 0.585270867049586 0.554142578592106 0.535605349852791];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS3D_p10(nTP)

%3D, alpha = 0.1

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 50 100 500];
slopeM = [0.001587591010958 0.000334860802492 0.000129492806027 0];
slopeP = [-0.001033885285809 -0.000409579718459 -0.000060967307118 0];
interseptM = [0.258642983727984 0.333806796235989 0.374880395529001 0.439626798542322];
interseptP = [0.644618927677761 0.613403649310276 0.57854240817611 0.548058754617302];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS3D_p05(nTP)

%3D, alpha = 0.05

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 50 80 500];
slopeM = [0.001897773535643 0.000390310143108 0.00014504631569 0];
slopeP = [-0.001281266825658 -0.000815336722161 -0.000079326637385 0];
interseptM = [0.217433363211081 0.307881166763189 0.356933932246803 0.429457090091877];
interseptP = [0.683369058696754 0.660072553521921 0.601191746739831 0.56152842804728];

% threshold curve evaluation

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

function [mssThreshNeg,mssThreshPos] = threshMSS3D_p01(nTP)

%3D, alpha = 0.01

%threshold curve parameters
turnPointsM = [20 50 200 500];
turnPointsP = [20 90 200 500];
slopeM = [0.00330864748993062 0.000607571513122873 0.000183704012624133 0];
slopeP = [-0.0011001388082068 -0.000401685416953439 -7.31736162754304e-05 0];
interseptM = [0.100955106735165 0.236008905575552 0.3207824056753 0.412634411987367];
interseptP = [0.747201850658008 0.684341045445206 0.618638685309604 0.582051877171889];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

%% threshold curve evaluation subfunction

function [mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,...
    turnPointsP,slopeM,slopeP,interseptM,interseptP,nTP)

%confined diffusion threshold
mssThreshNeg = NaN(1,turnPointsM(1)-1);
for i = 1 : length(turnPointsM)-1
    x = turnPointsM(i) : turnPointsM(i+1)-1;
    mssThreshNeg = [mssThreshNeg slopeM(i)*x+interseptM(i)]; %#ok<AGROW>
end
x = turnPointsM(end) : nTP;
mssThreshNeg = [mssThreshNeg slopeM(end)*x+interseptM(end)];

%directed diffusion threshold
mssThreshPos = NaN(1,turnPointsP(1)-1);
for i = 1 : length(turnPointsP)-1
    x = turnPointsP(i) : turnPointsP(i+1)-1;
    mssThreshPos = [mssThreshPos slopeP(i)*x+interseptP(i)]; %#ok<AGROW>
end
x = turnPointsP(end) : nTP;
mssThreshPos = [mssThreshPos slopeP(end)*x+interseptP(end)];

