function [ movie,metaData ] = readOME(fileAndPath)
%READOME read an OME tiff
% EHarry Jan 2012

if nargin < 1
    [file,dir]=uigetfile('*.ome.tiff'); % open file
    data = bfopenEdit(fullfile(dir,file)); % read file
else
    data = bfopenEdit(fileAndPath); % read file
    %[~,file] = fileparts(fileAndPath);
end

metaData = readMetaData(data{4});

% movie matrix goes [x,y,z,t,c]
movie = zeros(metaData.sizeX,metaData.sizeY,metaData.sizeZ,metaData.sizeT,metaData.sizeC);

if strcmp(metaData.dimOrder,'XYZCT') % only this dimORder for now
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

%%%% this is the old volocity metadata reader
%     function metaData = readMetaData(metadata)
%         metaString = metadata.get('<Pixels ID');
%
%         sizeX = regexp(metaString,'SizeX="\d+"','match');
%         sizeX = sizeX{1}(1:end-1);
%         metaData.sizeX = str2double(sizeX(8:end));
%
%         sizeY = regexp(metaString,'SizeY="\d+"','match');
%         sizeY = sizeY{1}(1:end-1);
%         metaData.sizeY = str2double(sizeY(8:end));
%
%         sizeZ = regexp(metaString,'SizeZ="\d+"','match');
%         sizeZ = sizeZ{1}(1:end-1);
%         metaData.sizeZ = str2double(sizeZ(8:end));
%
%         sizeC = regexp(metaString,'SizeC="\d+"','match');
%         sizeC = sizeC{1}(1:end-1);
%         metaData.sizeC = str2double(sizeC(8:end));
%
%         sizeT = regexp(metaString,'SizeT="\d+"','match');
%         sizeT = sizeT{1}(1:end-1);
%         metaData.sizeT = str2double(sizeT(8:end));
%
%         pixelXY = regexp(metaString,'PhysicalSizeX="[\d]+[.]?[\d]+"','match');
%         pixelXY = pixelXY{1}(1:end-1);
%         metaData.pixelXY = str2double(pixelXY(16:end));
%
%         pixelZ = regexp(metaString,'PhysicalSizeZ="[\d]+[.]?[\d]+"','match');
%         pixelZ = pixelZ{1}(1:end-1);
%         metaData.pixelZ = str2double(pixelZ(16:end));
%
%         dimOrder = regexp(metaString,'DimensionOrder="\w+"','match');
%         dimOrder = dimOrder{1}(1:end-1);
%         metaData.dimOrder = dimOrder(17:end);
%
%         metaString = metadata.get('<PlaneTiming DeltaT');
%
%         totalTime = regexp(metaString,'"[\d]+[.]?[\d]+"\sExp','match');
%         totalTime = totalTime{1}(1:end-5);
%         metaData.totalTime = str2double(totalTime(2:end));
%     end


%     function metaData = readMetaData(metadata,fileName) % this is for ome tifs written by imaris
%         metaString = metadata.get('<Pixels DimensionOrder');
%
%         sizeX = regexp(metaString,'SizeX = "\d+"','match');
%         sizeX = sizeX{1}(1:end-1);
%         metaData.sizeX = str2double(sizeX(10:end));
%
%         sizeY = regexp(metaString,'SizeY = "\d+"','match');
%         sizeY = sizeY{1}(1:end-1);
%         metaData.sizeY = str2double(sizeY(10:end));
%
%         sizeZ = regexp(metaString,'SizeZ = "\d+"','match');
%         sizeZ = sizeZ{1}(1:end-1);
%         metaData.sizeZ = str2double(sizeZ(10:end));
%
%         sizeC = regexp(metaString,'SizeC = "\d+"','match');
%         sizeC = sizeC{1}(1:end-1);
%         metaData.sizeC = str2double(sizeC(10:end));
%
%         sizeT = regexp(metaString,'SizeT = "\d+"','match');
%         sizeT = sizeT{1}(1:end-1);
%         metaData.sizeT = str2double(sizeT(10:end));
%
%         metaString = metadata.get('<Image Name');
%
%         pixelXY = regexp(metaString,'PixelSizeX = "[\d]+[.]?[\d]+"','match');
%         pixelXY = pixelXY{1}(1:end-1);
%         metaData.pixelXY = str2double(pixelXY(15:end));
%
%         pixelZ = regexp(metaString,'PixelSizeZ = "[\d]+[.]?[\d]+"','match');
%         pixelZ = pixelZ{1}(1:end-1);
%         metaData.pixelZ = str2double(pixelZ(15:end));
%
%         %         dimOrder = regexp(metaString,'DimensionOrder="\w+"','match');
%         %         dimOrder = dimOrder{1}(1:end-1);
%         metaData.dimOrder = 'XYZCT'; % only support this order for now
%
%         inD = inputdlg({['time lag for movie   "' fileName '"   in seconds']}); % have to ask for time lag
%
%         if isempty(inD)
%             error(['no input for time lag of movie   "' fileName '"   !']);
%         end
%
%         timeLag = str2double(inD{1});
%
%         if isnan(timeLag)
%             error(['not a valid time lag for movie   "' fileName '"   !']);
%         end
%
%         metaData.timeLag = timeLag;
%
%         metaString = metadata.get('<ChannelInfo ID');
%
%         wvl = regexp(metaString,'EmWave="\d+"','match');
%         wvl = wvl{1}(1:end-1);
%         metaData.wvl = str2double(wvl(9:end))/1000;
%     end

    function metaData = readMetaData(metadata) % this is for ome tifs written by volocity 6
        
        metaData.sizeX = str2double(metadata.getImage(0).getPixels.getSizeX.toString);
        metaData.sizeY = str2double(metadata.getImage(0).getPixels.getSizeY.toString);
        metaData.sizeZ = str2double(metadata.getImage(0).getPixels.getSizeZ.toString);
        metaData.sizeC = str2double(metadata.getImage(0).getPixels.getSizeC.toString);
        metaData.sizeT = str2double(metadata.getImage(0).getPixels.getSizeT.toString);
        
        metaData.pixelXY = str2double(metadata.getImage(0).getPixels.getPhysicalSizeX.toString);
        metaData.pixelZ = str2double(metadata.getImage(0).getPixels.getPhysicalSizeZ.toString);
        
        metaData.dimOrder = char(metadata.getImage(0).getPixels.getDimensionOrder);
        
        metaData.time = zeros(metaData.sizeZ.*metaData.sizeC.*metaData.sizeT,1);
        for i = 1:metaData.sizeZ.*metaData.sizeC.*metaData.sizeT
            metaData.time(i) = double(metadata.getImage(0).getPixels.getPlane(i-1).getDeltaT);
        end
        
        metaData.wvl = 0.5250; % only support gfp for now
    end

end

