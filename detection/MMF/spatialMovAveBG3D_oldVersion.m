function [bgMean,bgStd] = spatialMovAveBG3D_oldVersion(imageLast5,imageSizeX,imageSizeY,imageSizeZ)
% edit of spatialMovAveBG for 3d frames
% EHarry March 2012

%% ORINGIAL HEADER
% % %the function in its current form assigns blocks of 11x11 pixels the
% % %same background values, for the sake of speed

%define pixel limits where moving average can be calculated
startPixelX = 16;
endPixelX = max(imageSizeX - 15,startPixelX);
startPixelY = 16;
endPixelY = max(imageSizeY - 15,startPixelY);
startPixelZ = 16;
endPixelZ = max(imageSizeZ - 15,startPixelZ);

%allocate memory for output
bgMean = NaN(imageSizeX,imageSizeY,imageSizeZ);
bgStd = bgMean;

%go over all pixels within limits
for iPixelX = startPixelX : 11 : endPixelX
    for iPixelY = startPixelY : 11 : endPixelY
        for iPixelZ = startPixelZ : 11 : endPixelZ
            %get local image
            imageLocal = imageLast5(iPixelX-15:min(iPixelX+15,imageSizeX),iPixelY-15:min(iPixelY+15,imageSizeY),iPixelZ-15:min(iPixelZ+15,imageSizeZ),:);
            
            %estimate robust mean and std
            %first remove NaNs representing cropped regions
            imageLocal = imageLocal(~isnan(imageLocal));
            if ~isempty(imageLocal)
                [bgMean1,bgStd1] = robustMean(imageLocal(:));
                bgStd1 = max(bgStd1,eps);
            else
                bgMean1 = NaN;
                bgStd1 = NaN;
            end
            
            %put values in matrix representing image
            bgMean(iPixelX-5:iPixelX+5,iPixelY-5:iPixelY+5,iPixelZ-5:iPixelZ+5) = bgMean1;
            bgStd(iPixelX-5:iPixelX+5,iPixelY-5:iPixelY+5,iPixelZ-5:iPixelZ+5) = bgStd1;
        end
    end
end

%find limits of actual pixels filled up above
% firstFullX = find(~isnan(bgMean(:,startPixelY)),1,'first');
% lastFullX = find(~isnan(bgMean(:,startPixelY)),1,'last');
% firstFullY = find(~isnan(bgMean(startPixelX,:)),1,'first');
% lastFullY = find(~isnan(bgMean(startPixelX,:)),1,'last');
firstFullX = startPixelX - 5;
lastFullX = iPixelX + 5;
firstFullY = startPixelY - 5;
lastFullY = iPixelY + 5;
firstFullZ = startPixelZ - 5;
lastFullZ = iPixelZ + 5;

%patch the rest
% for iPixelZ = firstFullZ : lastFullZ
%     bgMean(1:firstFullX-1,1:firstFullY-1,iPixelZ) = bgMean(firstFullX,firstFullY,iPixelZ);
%     bgMean(lastFullX+1:end,lastFullY+1:end,iPixelZ) = bgMean(lastFullX,lastFullY,iPixelZ);
%     bgMean(1:firstFullX-1,lastFullY+1:end,iPixelZ) = bgMean(firstFullX,lastFullY,iPixelZ);
%     bgMean(lastFullX+1:end,1:firstFullY-1,iPixelZ) = bgMean(lastFullX,firstFullY,iPixelZ);
%
%     bgStd(1:firstFullX-1,1:firstFullY-1,iPixelZ) = bgStd(firstFullX,firstFullY,iPixelZ);
%     bgStd(lastFullX+1:end,lastFullY+1:end,iPixelZ) = bgStd(lastFullX,lastFullY,iPixelZ);
%     bgStd(1:firstFullX-1,lastFullY+1:end,iPixelZ) = bgStd(firstFullX,lastFullY,iPixelZ);
%     bgStd(lastFullX+1:end,1:firstFullY-1,iPixelZ) = bgStd(lastFullX,firstFullY,iPixelZ);
% end
%
% for iPixelY = firstFullY : lastFullY
%     bgMean(1:firstFullX-1,iPixelY,1:firstFullZ-1) = bgMean(firstFullX,iPixelY,firstFullZ);
%     bgMean(lastFullX+1:end,iPixelY,lastFullZ+1:end) = bgMean(lastFullX,iPixelY,lastFullZ);
%     bgMean(1:firstFullX-1,iPixelY,lastFullZ+1:end) = bgMean(firstFullX,iPixelY,lastFullZ);
%     bgMean(lastFullX+1:end,iPixelY,1:firstFullZ-1) = bgMean(lastFullX,iPixelY,firstFullZ);
%
%     bgStd(1:firstFullX-1,iPixelY,1:firstFullZ-1) = bgStd(firstFullX,iPixelY,firstFullZ);
%     bgStd(lastFullX+1:end,iPixelY,lastFullZ+1:end) = bgStd(lastFullX,iPixelY,lastFullZ);
%     bgStd(1:firstFullX-1,iPixelY,lastFullZ+1:end) = bgStd(firstFullX,iPixelY,lastFullZ);
%     bgStd(lastFullX+1:end,iPixelY,1:firstFullZ-1) = bgStd(lastFullX,iPixelY,firstFullZ);
% end
%
% for iPixelX = 1 : imageSizeX
%     bgMean(iPixelX,1:firstFullY-1,1:firstFullZ-1) = bgMean(iPixelX,firstFullY,firstFullZ);
%     bgMean(iPixelX,lastFullY+1:end,lastFullZ+1:end) = bgMean(iPixelX,lastFullY,lastFullZ);
%     bgMean(iPixelX,1:firstFullY-1,lastFullZ+1:end) = bgMean(iPixelX,firstFullY,lastFullZ);
%     bgMean(iPixelX,lastFullY+1:end,1:firstFullZ-1) = bgMean(iPixelX,lastFullY,firstFullZ);
%
%     bgStd(iPixelX,1:firstFullY-1,1:firstFullZ-1) = bgStd(iPixelX,firstFullY,firstFullZ);
%     bgStd(iPixelX,lastFullY+1:end,lastFullZ+1:end) = bgStd(iPixelX,lastFullY,lastFullZ);
%     bgStd(iPixelX,1:firstFullY-1,lastFullZ+1:end) = bgStd(iPixelX,firstFullY,lastFullZ);
%     bgStd(iPixelX,lastFullY+1:end,1:firstFullZ-1) = bgStd(iPixelX,lastFullY,firstFullZ);
% end

for iPixelZ = firstFullZ : lastFullZ
    for iPixelY = firstFullY : lastFullY
        bgMean(1:firstFullX-1,iPixelY,1:iPixelZ) = bgMean(firstFullX,iPixelY,iPixelZ);
        bgMean(lastFullX+1:end,iPixelY,iPixelZ) = bgMean(lastFullX,iPixelY,iPixelZ);
        bgMean(1:firstFullX-1,iPixelY,iPixelZ) = bgMean(firstFullX,iPixelY,iPixelZ);
        bgMean(lastFullX+1:end,iPixelY,iPixelZ) = bgMean(lastFullX,iPixelY,iPixelZ);
        
        bgStd(1:firstFullX-1,iPixelY,iPixelZ) = bgStd(firstFullX,iPixelY,iPixelZ);
        bgStd(lastFullX+1:end,iPixelY,iPixelZ) = bgStd(lastFullX,iPixelY,iPixelZ);
        bgStd(1:firstFullX-1,iPixelY,iPixelZ) = bgStd(firstFullX,iPixelY,iPixelZ);
        bgStd(lastFullX+1:end,iPixelY,iPixelZ) = bgStd(lastFullX,iPixelY,iPixelZ);
    end
    
    for iPixelX = 1 : imageSizeX
        bgMean(iPixelX,1:firstFullY-1,iPixelZ) = bgMean(iPixelX,firstFullY,iPixelZ);
        bgMean(iPixelX,lastFullY+1:end,iPixelZ) = bgMean(iPixelX,lastFullY,iPixelZ);
        bgMean(iPixelX,1:firstFullY-1,iPixelZ) = bgMean(iPixelX,firstFullY,iPixelZ);
        bgMean(iPixelX,lastFullY+1:end,iPixelZ) = bgMean(iPixelX,lastFullY,iPixelZ);
        
        bgStd(iPixelX,1:firstFullY-1,iPixelZ) = bgStd(iPixelX,firstFullY,iPixelZ);
        bgStd(iPixelX,lastFullY+1:end,iPixelZ) = bgStd(iPixelX,lastFullY,iPixelZ);
        bgStd(iPixelX,1:firstFullY-1,iPixelZ) = bgStd(iPixelX,firstFullY,iPixelZ);
        bgStd(iPixelX,lastFullY+1:end,iPixelZ) = bgStd(iPixelX,lastFullY,iPixelZ);
    end
end

for iPixelY = 1 : imageSizeY
    for iPixelX = 1: imageSizeX
        bgMean(iPixelX,iPixelY,1:firstFullZ-1) = bgMean(iPixelX,iPixelY,firstFullZ);
        bgMean(iPixelX,iPixelY,lastFullZ+1:end) = bgMean(iPixelX,iPixelY,lastFullZ);
        bgMean(iPixelX,iPixelY,lastFullZ+1:end) = bgMean(iPixelX,iPixelY,lastFullZ);
        bgMean(iPixelX,iPixelY,1:firstFullZ-1) = bgMean(iPixelX,iPixelY,firstFullZ);
        
        bgStd(iPixelX,iPixelY,1:firstFullZ-1) = bgStd(iPixelX,iPixelY,firstFullZ);
        bgStd(iPixelX,iPixelY,lastFullZ+1:end) = bgStd(iPixelX,iPixelY,lastFullZ);
        bgStd(iPixelX,iPixelY,lastFullZ+1:end) = bgStd(iPixelX,iPixelY,lastFullZ);
        bgStd(iPixelX,iPixelY,1:firstFullZ-1) = bgStd(iPixelX,iPixelY,firstFullZ);
    end
end

% for iPixelX = 1 : imageSizeX
%     bgMean(iPixelX,1:firstFullY-1,1:firstFullZ-1) = bgMean(iPixelX,firstFullY,firstFullZ);
%     bgMean(iPixelX,lastFullY+1:end,lastFullZ+1:end) = bgMean(iPixelX,lastFullY,lastFullZ);
%     bgMean(iPixelX,1:firstFullY-1,lastFullZ+1:end) = bgMean(iPixelX,firstFullY,lastFullZ);
%     bgMean(iPixelX,lastFullY+1:end,1:firstFullZ-1) = bgMean(iPixelX,lastFullY,firstFullZ);
%
%     bgStd(iPixelX,1:firstFullY-1,1:firstFullZ-1) = bgStd(iPixelX,firstFullY,firstFullZ);
%     bgStd(iPixelX,lastFullY+1:end,lastFullZ+1:end) = bgStd(iPixelX,lastFullY,lastFullZ);
%     bgStd(iPixelX,1:firstFullY-1,lastFullZ+1:end) = bgStd(iPixelX,firstFullY,lastFullZ);
%     bgStd(iPixelX,lastFullY+1:end,1:firstFullZ-1) = bgStd(iPixelX,lastFullY,firstFullZ);
% end

