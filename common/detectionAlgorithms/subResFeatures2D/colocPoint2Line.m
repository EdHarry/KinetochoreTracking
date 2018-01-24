function [fracPointsOnLines,pValue,numPoints,fracPointsOnLinesRand] = ...
    colocPoint2Line(lineImage,pointImage,lineDetectInput,pointDetectInput,...
    numRepTestRand,doPlot)

%% input

%line detection parameters
scales = lineDetectInput.scales;
lineType = lineDetectInput.lineType;
conf = lineDetectInput.conf;
dilateW = lineDetectInput.dilateW;

%point detection parameters
alphaLocMax = pointDetectInput.alphaLocMax;

%get image size
[imgSizeX,imgSizeY] = size(pointImage);

%% line detection

%detect lines in lineImage
img.data = lineImage;
img.perm = 'M';
respRaw = imLineDetect(img,scales,lineType,conf);

%keep only significant lines - HOW???
% respVal = respRaw(respRaw~=0);
% respThresh = prctile(respVal,5);
respThresh = 0;
resp = respRaw;
resp(resp<=respThresh) = 0;
resp(resp>respThresh) = 1;

%dilate resp by requested width
SE = strel('square',dilateW);
respDilate = imdilate(resp,SE);

%% point detection

%filter pointImage
pointImageF = Gauss2D(pointImage,1);

%get background information
[bgMean,bgStd] = spatialMovAveBG(pointImageF,size(pointImageF,1),size(pointImageF,2));

%find local maxima in filtered image
fImg = locmax2d(pointImageF,[3 3],1);

%get positions and amplitudes of local maxima
[localMaxPosX,localMaxPosY,localMaxAmp] = find(fImg);
localMax1DIndx = find(fImg(:));

%get background values corresponding to local maxima
bgMeanMax = bgMean(localMax1DIndx);
bgStdMax = bgStd(localMax1DIndx);

%calculate the p-value corresponding to the local maxima's amplitudes
%assume that background intensity in filtered image is normally
%distributed with mean bgMeanMax and standard deviation bgStdMax
pValue = 1 - normcdf(localMaxAmp,bgMeanMax,bgStdMax);

%retain only those maxima with significant amplitude
keepMax = find(pValue < alphaLocMax);
localMaxPosX = localMaxPosX(keepMax);
localMaxPosY = localMaxPosY(keepMax);
positions = [localMaxPosX localMaxPosY];

%convert positions from a 2D index to a 1D index
positions1D = (positions(:,2)-1)*imgSizeX + positions(:,1);

%get total number of detected points
numPoints = size(positions,1);

%% colocalization

%get number of points that colocalize with lines
fracPointsOnLines = length(find(respDilate(positions1D)==1))/numPoints;

%% test colocalization with a random distribution

%get how many times to repeat calculation
numRep = numRepTestRand;

%generate numPoints with random positions, numRep times
randomPos = ceil(rand(numPoints,numRep)*imgSizeX*imgSizeY);

%get number of randomly distributed points that colocalize with lines
fracPointsOnLinesRand = NaN(numRep,1);
for iRep = 1 : numRep
    fracPointsOnLinesRand(iRep) = length(find(respDilate(randomPos...
        (:,iRep))==1))/numPoints;
end

%% assuming a normal distribution, get observed value's p-value

meanRand = mean(fracPointsOnLinesRand);
stdRand = std(fracPointsOnLinesRand);

if fracPointsOnLines > meanRand
    pValue = 1 - normcdf(fracPointsOnLines,meanRand,stdRand);
else
    pValue = normcdf(fracPointsOnLines,meanRand,stdRand);
end

%% plot if requested
if doPlot

    %plot detected points on top of pointImage
    plotImageWithFeatures(pointImage,[positions(:,2) positions(:,1)]);
    
    %plot detected lines on top of lineImage
    [linePosX,linePosY] = find(resp~=0);
    plotImageWithFeatures(lineImage,[linePosY linePosX]);
    
    %plot detected points on top of detected lines
    plotImageWithFeatures(respDilate,[positions(:,2) positions(:,1)]);
    
end

%% ~~~ the end ~~~


