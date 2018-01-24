function headerOME = convertOMEMetaDataToMakiHeaderVolocityVersion( metaData )
%CONVERTOMEMETADATATOMAKIHEADER convert OME metadata to maki compatable
%header
% EHarry October 2011

headerOME.pixelX = metaData.pixelXY;
headerOME.pixelY = metaData.pixelXY;
headerOME.pixelZ = metaData.pixelZ;
headerOME.numCols = metaData.sizeY;
headerOME.numRows = metaData.sizeX;
headerOME.numZSlices = metaData.sizeZ;
headerOME.numTimepoints = metaData.sizeT;
headerOME.numWvs = 1;
headerOME.zwtOrder = 'ztw';
headerOME.Time = zeros(metaData.sizeZ.*metaData.sizeT,1);
headerOME.Time(1:end) = deal(metaData.totalTime./(metaData.sizeT-1));
headerOME.timestamp = zeros(metaData.sizeZ,metaData.sizeT);
headerOME.timestamp(1:end) = deal(metaData.totalTime./(metaData.sizeT-1));
headerOME.lensID = 10007;
headerOME.wvl = 0.5250;

end

