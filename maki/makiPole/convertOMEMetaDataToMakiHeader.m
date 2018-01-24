function headerOME = convertOMEMetaDataToMakiHeader( metaData )
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
headerOME.numWvs = metaData.sizeC;
headerOME.zwtOrder = 'ztw';
% headerOME.timestamp = zeros(metaData.sizeZ,metaData.sizeT);
% for i = 0:metaData.sizeT-1
%     headerOME.timestamp(:,i+1) = i.*metaData.timeLag.*ones(metaData.sizeZ,1);
% end
% headerOME.Time = reshape(headerOME.timestamp,metaData.sizeZ.*metaData.sizeT,1);
headerOME.Time = metaData.time;
timeStamp = reshape(metaData.time,metaData.sizeZ,metaData.sizeC,metaData.sizeT);
headerOME.timestamp = permute(timeStamp,[1 3 2]);
headerOME.lensID = 10007; % assume this lens ID (100x) for now.
headerOME.wvl = metaData.wvl;
end

