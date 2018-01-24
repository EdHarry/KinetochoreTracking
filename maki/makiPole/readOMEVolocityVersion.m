function [ movie,metaData ] = readOMEVolocityVersion(fileAndPath)
% EHarry Feb 2012

if nargin < 1
    [file,dir]=uigetfile('*'); % open file
    data = bfopenEdit(fullfile(dir,file)); % read file
else
    data = bfopenEdit(fileAndPath); % read file
end

metaData = readMetaData(data{2});

% movie matrix goes [x,y,z,t,c]
movie = zeros(metaData.sizeX,metaData.sizeY,metaData.sizeZ,metaData.sizeT,metaData.sizeC);

if strcmp(metaData.dimOrder,'XYZCT') % only this dimORder for now since this is the format of volocity ome tiffs
    count=1;
    for t = 1:metaData.sizeT
        for c = 1:metaData.sizeC
            for z = 1:metaData.sizeZ
                movie(:,:,z,t,c) = data{1}{count,1}'; % data goes (y,x) hence the transpose
                count=count+1;
            end
        end
    end
end


%% SUBFUNCTIONS

    function metaData = readMetaData(metadata)
        metaString = metadata.get('<Pixels ID');
        
        sizeX = regexp(metaString,'SizeX="\d+"','match');
        sizeX = sizeX{1}(1:end-1);
        metaData.sizeX = str2double(sizeX(8:end));
        
        sizeY = regexp(metaString,'SizeY="\d+"','match');
        sizeY = sizeY{1}(1:end-1);
        metaData.sizeY = str2double(sizeY(8:end));
        
        sizeZ = regexp(metaString,'SizeZ="\d+"','match');
        sizeZ = sizeZ{1}(1:end-1);
        metaData.sizeZ = str2double(sizeZ(8:end));
        
        sizeC = regexp(metaString,'SizeC="\d+"','match');
        sizeC = sizeC{1}(1:end-1);
        metaData.sizeC = str2double(sizeC(8:end));
        
        sizeT = regexp(metaString,'SizeT="\d+"','match');
        sizeT = sizeT{1}(1:end-1);
        metaData.sizeT = str2double(sizeT(8:end));
        
        pixelXY = regexp(metaString,'PhysicalSizeX="[\d]+[.]?[\d]+"','match');
        pixelXY = pixelXY{1}(1:end-1);
        metaData.pixelXY = str2double(pixelXY(16:end));
        
        pixelZ = regexp(metaString,'PhysicalSizeZ="[\d]+[.]?[\d]+"','match');
        pixelZ = pixelZ{1}(1:end-1);
        metaData.pixelZ = str2double(pixelZ(16:end));
        
        dimOrder = regexp(metaString,'DimensionOrder="\w+"','match');
        dimOrder = dimOrder{1}(1:end-1);
        metaData.dimOrder = dimOrder(17:end);
        
        metaString = metadata.get('<PlaneTiming DeltaT');
        
        totalTime = regexp(metaString,'"[\d]+[.]?[\d]+"\sExp','match');
        totalTime = totalTime{1}(1:end-5);
        metaData.totalTime = str2double(totalTime(2:end));
    end

end

