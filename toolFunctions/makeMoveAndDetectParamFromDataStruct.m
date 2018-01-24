function [ movieParam , detectParam ] = makeMoveAndDetectParamFromDataStruct( dataStruct, trackingChNo , decon , crop , firstFrame , lastFrame, filterPrm, integWindow, cutOffMethod, testAlpha, alphaLocMax )
%MAKEMOVEANDDETECTPARAMFROMDATASTRUCT makes a movieParam and detectParam
%struct for MMF from a maki dataStruct
%
%   dataStruct:   maki dataStruct
%
%
%   OPTIONAL (can either not be inputed or inputed as []):
%
%   IMAGE PARAMETERS:
%   trackingChNo : index to the channel to be tracked
%   decon :        use deconvolved version of the movie or not (0/1)
%   crop:          crop matrix for the movie (maki style)
%   firstFrame:    index of the first frame to use
%   lastFrame:     index of the last frame to use
%
%                  image parameters will be taken from the dataStruct if
%                  not inputed (except trackingChNo, not in dataStruct yet, defaults to 1 for now and firstFrame, defaults to 1)
%
%
%   DETECTION PARMETERS:
%   filterPrm:      filtering parameters, maki style, defaults to the ones
%                   in the dataStruct 
%   integWindow:    localMaxima intergration windows, defaults to 0
%   
%   cutOffMethod:   localMax cutoff method, defaults to 'uniModal'
%   testAlpha:      MMF alpha value cutoffs, defaults to struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0)
%   alphaLocMax:    localMax alpha cutoff (only relavent if cutOffMethod is
%                   'pValue'), defaults to 1e-6 
%
%
%   EHarry March 2012


%% check input
if nargin < 1 || isempty(dataStruct)
    error('dataStruct must be input');
end

if nargin < 2
    trackingChNo = [];
end

if nargin < 3
    decon = [];
end

if nargin < 4
    crop = [];
end

if nargin < 5
    firstFrame = [];
end

if nargin < 6
    lastFrame = [];
end

if nargin < 7
    filterPrm = [];
end

if nargin < 8
    integWindow = [];
end

if nargin < 9
    cutOffMethod = [];
end

if nargin < 10
    testAlpha = [];
end

if nargin < 11
    alphaLocMax = [];
end

%% define defaults, some from the dataStruct

% image parameters
trackingChNo_def = 1;

if isfield(dataStruct.dataProperties,'decon')
    decon_def_tmp = dataStruct.dataProperties.decon;
else
    decon_def_tmp = 0;
end
decon_def = decon_def_tmp;

if isfield(dataStruct.dataProperties,'crop')
    crop_def_tmp = dataStruct.dataProperties.crop;
else
    crop_def_tmp = [];
end
crop_def = crop_def_tmp;

firstFrame_def = 1;

lastFrame_def = dataStruct.dataProperties.movieSize(4);

% detection param
filterPrm_def = dataStruct.dataProperties.FILTERPRM;

integWindow_def = 0;

cutOffMethod_def = 'uniModal';

testAlpha_def = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0);

alphaLocMax_def = 1e-6;





%% moveParam
movieParam.imageName = fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName);

if isempty(trackingChNo)
    movieParam.imageCh = trackingChNo_def;
else
    movieParam.imageCh = trackingChNo;
end

if isempty(decon)
    movieParam.imageDecon = decon_def;
else
    movieParam.imageDecon = decon;
end

if isempty(crop)
    movieParam.imageCrop = crop_def;
else
    movieParam.imageCrop = crop;
end

if isempty(firstFrame)
    movieParam.firstImageNum = firstFrame_def;
else
    movieParam.firstImageNum = firstFrame;
end

if isempty(lastFrame)
    movieParam.lastImageNum = lastFrame_def;
else
    movieParam.lastImageNum = lastFrame;
end


%% detectParam
if isempty(filterPrm)
    detectParam.filterPrm = filterPrm_def;
else
    detectParam.filterPrm = filterPrm;
end

if isempty(integWindow)
    detectParam.integWindow = integWindow_def;
else
    detectParam.integWindow = integWindow;
end

if isempty(cutOffMethod)
    detectParam.cutOffMethod = cutOffMethod_def;
else
    detectParam.cutOffMethod = cutOffMethod;
end

if isempty(testAlpha)
    detectParam.testAlpha = testAlpha_def;
else
    detectParam.testAlpha = testAlpha;
end

if isempty(alphaLocMax)
    detectParam.alphaLocMax = alphaLocMax_def;
else
    detectParam.alphaLocMax = alphaLocMax;
end

end

