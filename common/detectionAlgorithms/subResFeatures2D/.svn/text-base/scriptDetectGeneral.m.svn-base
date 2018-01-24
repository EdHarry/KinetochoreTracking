
%% movie information
movieParam.imageDir = '/mnt/sickkids/Yoav/2009_05_06_monodisperse_QDots_blinking_calibration/8ms/150sens/cell_02/images/'; %directory where images are
movieParam.filenameBase = 'monodisperse_new_camera_calibration_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 100; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 2.1; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',5e-5,'alphaA',5e-6,'alphaD',5e-6,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.01; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

%% save results
saveResults.dir = '/mnt/sickkids/Yoav/2009_05_06_monodisperse_QDots_blinking_calibration/8ms/150sens/cell_02/analysis/'; %directory where to save input and output
saveResults.filename = 'detectionTest3.mat'; %name of file where input and output are saved
% saveResults = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
