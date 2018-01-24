function dataStruct = makiGaussianFitting_mlSparse( dataStruct)
%MAKIGAUSSIANFITTING get initCoords from gaussian fitting
%   EHarry May 2012


% first make detectionParams
[ movieParam , detectParam ] = makeMoveAndDetectParamFromDataStruct( dataStruct, 1 , [] , [] , [] , [], [], 0, [] , struct('alphaR',0.05,'alphaA',0.01,'alphaD',0.01,'alphaF',0), 1e-7 );


% then run fitting
[movieInfo,exceptions,localMaxima,~,psfSigma] = detectSpots_newVersion_mlSparse_reCluster(movieParam,detectParam,'no');

% convert into initCoords
dataStruct.initCoord = convertMovieInfo2InitCoord( movieInfo,dataStruct.dataProperties );

% save extra stuff
dataStruct.initCoord(1).exceptions = exceptions;
dataStruct.initCoord(1).localMaxima = localMaxima;
dataStruct.dataProperties.psfSigma = psfSigma;

% initCoord depends on dataProperites
dependencies = struct('dataProperties',[]);

dataPropName = dataStruct.dataPropertiesName;
dataPropV = getVersion(dataPropName);

dependencies.dataProperties = dataPropV;

% save into the first initCoords
dataStruct.initCoord(1).dependencies = dependencies;
end

