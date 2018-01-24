function [bgMean,bgStd] = spatialMovAveBG3D_newVersion(imageLast5,imageSizeX,imageSizeY,imageSizeZ,verbose)
% edit of spatialMovAveBG for 3d frames
% EHarry March 2012

%% ORINGIAL HEADER
% % %the function in its current form assigns blocks of 11x11x11 pixels the
% % %same background values, for the sake of speed

%% OLD CODE
% % %define pixel limits where moving average can be calculated
% % startPixelX = 16;
% % endPixelX = max(imageSizeX - 15,startPixelX);
% % startPixelY = 16;
% % endPixelY = max(imageSizeY - 15,startPixelY);
% % startPixelZ = 16;
% % endPixelZ = max(imageSizeZ - 15,startPixelZ);
% %
% % %allocate memory for output
% % bgMean = NaN(imageSizeX,imageSizeY,imageSizeZ);
% % bgStd = bgMean;
% %
% % %go over all pixels within limits
% % for iPixelX = startPixelX : 11 : endPixelX
% %     for iPixelY = startPixelY : 11 : endPixelY
% %         for iPixelZ = startPixelZ : 11 : endPixelZ
% %             %get local image
% %             imageLocal = imageLast5(iPixelX-15:min(iPixelX+15,imageSizeX),iPixelY-15:min(iPixelY+15,imageSizeY),iPixelZ-15:min(iPixelZ+15,imageSizeZ),:);
% %
% %             %estimate robust mean and std
% %             %first remove NaNs representing cropped regions
% %             imageLocal = imageLocal(~isnan(imageLocal));
% %             if ~isempty(imageLocal)
% %                 [bgMean1,bgStd1] = robustMean(imageLocal(:));
% %                 bgStd1 = max(bgStd1,eps);
% %             else
% %                 bgMean1 = NaN;
% %                 bgStd1 = NaN;
% %             end
% %
% %             %put values in matrix representing image
% %             bgMean(iPixelX-5:iPixelX+5,iPixelY-5:iPixelY+5,iPixelZ-5:iPixelZ+5) = bgMean1;
% %             bgStd(iPixelX-5:iPixelX+5,iPixelY-5:iPixelY+5,iPixelZ-5:iPixelZ+5) = bgStd1;
% %         end
% %     end
% % end
% %
% % %find limits of actual pixels filled up above
% % % firstFullX = find(~isnan(bgMean(:,startPixelY)),1,'first');
% % % lastFullX = find(~isnan(bgMean(:,startPixelY)),1,'last');
% % % firstFullY = find(~isnan(bgMean(startPixelX,:)),1,'first');
% % % lastFullY = find(~isnan(bgMean(startPixelX,:)),1,'last');
% % firstFullX = startPixelX - 5;
% % lastFullX = iPixelX + 5;
% % firstFullY = startPixelY - 5;
% % lastFullY = iPixelY + 5;
% % firstFullZ = startPixelZ - 5;
% % lastFullZ = iPixelZ + 5;
% %
% % %patch the rest
% % % for iPixelZ = firstFullZ : lastFullZ
% % %     bgMean(1:firstFullX-1,1:firstFullY-1,iPixelZ) = bgMean(firstFullX,firstFullY,iPixelZ);
% % %     bgMean(lastFullX+1:end,lastFullY+1:end,iPixelZ) = bgMean(lastFullX,lastFullY,iPixelZ);
% % %     bgMean(1:firstFullX-1,lastFullY+1:end,iPixelZ) = bgMean(firstFullX,lastFullY,iPixelZ);
% % %     bgMean(lastFullX+1:end,1:firstFullY-1,iPixelZ) = bgMean(lastFullX,firstFullY,iPixelZ);
% % %
% % %     bgStd(1:firstFullX-1,1:firstFullY-1,iPixelZ) = bgStd(firstFullX,firstFullY,iPixelZ);
% % %     bgStd(lastFullX+1:end,lastFullY+1:end,iPixelZ) = bgStd(lastFullX,lastFullY,iPixelZ);
% % %     bgStd(1:firstFullX-1,lastFullY+1:end,iPixelZ) = bgStd(firstFullX,lastFullY,iPixelZ);
% % %     bgStd(lastFullX+1:end,1:firstFullY-1,iPixelZ) = bgStd(lastFullX,firstFullY,iPixelZ);
% % % end
% % %
% % % for iPixelY = firstFullY : lastFullY
% % %     bgMean(1:firstFullX-1,iPixelY,1:firstFullZ-1) = bgMean(firstFullX,iPixelY,firstFullZ);
% % %     bgMean(lastFullX+1:end,iPixelY,lastFullZ+1:end) = bgMean(lastFullX,iPixelY,lastFullZ);
% % %     bgMean(1:firstFullX-1,iPixelY,lastFullZ+1:end) = bgMean(firstFullX,iPixelY,lastFullZ);
% % %     bgMean(lastFullX+1:end,iPixelY,1:firstFullZ-1) = bgMean(lastFullX,iPixelY,firstFullZ);
% % %
% % %     bgStd(1:firstFullX-1,iPixelY,1:firstFullZ-1) = bgStd(firstFullX,iPixelY,firstFullZ);
% % %     bgStd(lastFullX+1:end,iPixelY,lastFullZ+1:end) = bgStd(lastFullX,iPixelY,lastFullZ);
% % %     bgStd(1:firstFullX-1,iPixelY,lastFullZ+1:end) = bgStd(firstFullX,iPixelY,lastFullZ);
% % %     bgStd(lastFullX+1:end,iPixelY,1:firstFullZ-1) = bgStd(lastFullX,iPixelY,firstFullZ);
% % % end
% % %
% % % for iPixelX = 1 : imageSizeX
% % %     bgMean(iPixelX,1:firstFullY-1,1:firstFullZ-1) = bgMean(iPixelX,firstFullY,firstFullZ);
% % %     bgMean(iPixelX,lastFullY+1:end,lastFullZ+1:end) = bgMean(iPixelX,lastFullY,lastFullZ);
% % %     bgMean(iPixelX,1:firstFullY-1,lastFullZ+1:end) = bgMean(iPixelX,firstFullY,lastFullZ);
% % %     bgMean(iPixelX,lastFullY+1:end,1:firstFullZ-1) = bgMean(iPixelX,lastFullY,firstFullZ);
% % %
% % %     bgStd(iPixelX,1:firstFullY-1,1:firstFullZ-1) = bgStd(iPixelX,firstFullY,firstFullZ);
% % %     bgStd(iPixelX,lastFullY+1:end,lastFullZ+1:end) = bgStd(iPixelX,lastFullY,lastFullZ);
% % %     bgStd(iPixelX,1:firstFullY-1,lastFullZ+1:end) = bgStd(iPixelX,firstFullY,lastFullZ);
% % %     bgStd(iPixelX,lastFullY+1:end,1:firstFullZ-1) = bgStd(iPixelX,lastFullY,firstFullZ);
% % % end
% %
% % for iPixelZ = firstFullZ : lastFullZ
% %     for iPixelY = firstFullY : lastFullY
% %         bgMean(1:firstFullX-1,iPixelY,1:iPixelZ) = bgMean(firstFullX,iPixelY,iPixelZ);
% %         bgMean(lastFullX+1:end,iPixelY,iPixelZ) = bgMean(lastFullX,iPixelY,iPixelZ);
% %         bgMean(1:firstFullX-1,iPixelY,iPixelZ) = bgMean(firstFullX,iPixelY,iPixelZ);
% %         bgMean(lastFullX+1:end,iPixelY,iPixelZ) = bgMean(lastFullX,iPixelY,iPixelZ);
% %
% %         bgStd(1:firstFullX-1,iPixelY,iPixelZ) = bgStd(firstFullX,iPixelY,iPixelZ);
% %         bgStd(lastFullX+1:end,iPixelY,iPixelZ) = bgStd(lastFullX,iPixelY,iPixelZ);
% %         bgStd(1:firstFullX-1,iPixelY,iPixelZ) = bgStd(firstFullX,iPixelY,iPixelZ);
% %         bgStd(lastFullX+1:end,iPixelY,iPixelZ) = bgStd(lastFullX,iPixelY,iPixelZ);
% %     end
% %
% %     for iPixelX = 1 : imageSizeX
% %         bgMean(iPixelX,1:firstFullY-1,iPixelZ) = bgMean(iPixelX,firstFullY,iPixelZ);
% %         bgMean(iPixelX,lastFullY+1:end,iPixelZ) = bgMean(iPixelX,lastFullY,iPixelZ);
% %         bgMean(iPixelX,1:firstFullY-1,iPixelZ) = bgMean(iPixelX,firstFullY,iPixelZ);
% %         bgMean(iPixelX,lastFullY+1:end,iPixelZ) = bgMean(iPixelX,lastFullY,iPixelZ);
% %
% %         bgStd(iPixelX,1:firstFullY-1,iPixelZ) = bgStd(iPixelX,firstFullY,iPixelZ);
% %         bgStd(iPixelX,lastFullY+1:end,iPixelZ) = bgStd(iPixelX,lastFullY,iPixelZ);
% %         bgStd(iPixelX,1:firstFullY-1,iPixelZ) = bgStd(iPixelX,firstFullY,iPixelZ);
% %         bgStd(iPixelX,lastFullY+1:end,iPixelZ) = bgStd(iPixelX,lastFullY,iPixelZ);
% %     end
% % end
% %
% % for iPixelY = 1 : imageSizeY
% %     for iPixelX = 1: imageSizeX
% %         bgMean(iPixelX,iPixelY,1:firstFullZ-1) = bgMean(iPixelX,iPixelY,firstFullZ);
% %         bgMean(iPixelX,iPixelY,lastFullZ+1:end) = bgMean(iPixelX,iPixelY,lastFullZ);
% %         bgMean(iPixelX,iPixelY,lastFullZ+1:end) = bgMean(iPixelX,iPixelY,lastFullZ);
% %         bgMean(iPixelX,iPixelY,1:firstFullZ-1) = bgMean(iPixelX,iPixelY,firstFullZ);
% %
% %         bgStd(iPixelX,iPixelY,1:firstFullZ-1) = bgStd(iPixelX,iPixelY,firstFullZ);
% %         bgStd(iPixelX,iPixelY,lastFullZ+1:end) = bgStd(iPixelX,iPixelY,lastFullZ);
% %         bgStd(iPixelX,iPixelY,lastFullZ+1:end) = bgStd(iPixelX,iPixelY,lastFullZ);
% %         bgStd(iPixelX,iPixelY,1:firstFullZ-1) = bgStd(iPixelX,iPixelY,firstFullZ);
% %     end
% % end
% %
% % % for iPixelX = 1 : imageSizeX
% % %     bgMean(iPixelX,1:firstFullY-1,1:firstFullZ-1) = bgMean(iPixelX,firstFullY,firstFullZ);
% % %     bgMean(iPixelX,lastFullY+1:end,lastFullZ+1:end) = bgMean(iPixelX,lastFullY,lastFullZ);
% % %     bgMean(iPixelX,1:firstFullY-1,lastFullZ+1:end) = bgMean(iPixelX,firstFullY,lastFullZ);
% % %     bgMean(iPixelX,lastFullY+1:end,1:firstFullZ-1) = bgMean(iPixelX,lastFullY,firstFullZ);
% % %
% % %     bgStd(iPixelX,1:firstFullY-1,1:firstFullZ-1) = bgStd(iPixelX,firstFullY,firstFullZ);
% % %     bgStd(iPixelX,lastFullY+1:end,lastFullZ+1:end) = bgStd(iPixelX,lastFullY,lastFullZ);
% % %     bgStd(iPixelX,1:firstFullY-1,lastFullZ+1:end) = bgStd(iPixelX,firstFullY,lastFullZ);
% % %     bgStd(iPixelX,lastFullY+1:end,1:firstFullZ-1) = bgStd(iPixelX,lastFullY,firstFullZ);
% end


if nargin < 5 || isempty(verbose)
    verbose = 0;
end

% average over 11x11x3 pixel cubes
% less in z because of the expected difference in pixel sizes
% these have to be odd
rangeX = 11;
rangeY = 11;
rangeZ = 3;


bgMean = NaN(imageSizeX,imageSizeY,imageSizeZ);
bgStd = bgMean;

if verbose
    progressText(0,'Calculating local background at {x,y,z}...');
end

% timeCount = 0;
% timeTotal = imageSizeX*imageSizeY*imageSizeZ;

% new, assign blocks of pixels of size (rangeX,rangeY,rangeZ)/2 to the same
% values for speed

loopRangeX = rangeX;
loopRangeY = rangeY;
loopRangeZ = rangeZ;

% loopRangeX = (rangeX-1)/2;
% loopRangeY = (rangeY-1)/2;
% loopRangeZ = (rangeZ-1)/2;

timeCount = 0;
timeTotal = ceil(imageSizeX/loopRangeX)*ceil(imageSizeY/loopRangeY)*ceil(imageSizeZ/loopRangeZ);


% for iX = 1:rangeX:imageSizeX
%     for iY = 1:rangeY:imageSizeY
%         for iZ = 1:rangeZ:imageSizeZ

for iX = 1:loopRangeX:imageSizeX
    for iY = 1:loopRangeY:imageSizeY
        for iZ = 1:loopRangeZ:imageSizeZ
            
            % get middle pixels of the range
            iXm = min(iX+((rangeX-1)/2),imageSizeX);
            iYm = min(iY+((rangeY-1)/2),imageSizeY);
            iZm = min(iZ+((rangeZ-1)/2),imageSizeZ);
            
            %             % get the limits for the box for the current pixel
            %             xMin = max(iX-((rangeX-1)/2),1);
            %             xMax = min(xMin+rangeX-1,imageSizeX);
            %             xMin = max(xMax-rangeX+1,1);
            %             yMin = max(iY-((rangeY-1)/2),1);
            %             yMax = min(yMin+rangeY-1,imageSizeY);
            %             yMin = max(yMax-rangeY+1,1);
            %             zMin = max(iZ-((rangeZ-1)/2),1);
            %             zMax = min(zMin+rangeZ-1,imageSizeZ);
            %             zMin = max(zMax-rangeZ+1,1);
            
            % get the limits for the box for the current pixel
            xMin = max(iXm-((rangeX-1)/2),1);
            xMax = min(xMin+rangeX-1,imageSizeX);
            xMin = max(xMax-rangeX+1,1);
            yMin = max(iYm-((rangeY-1)/2),1);
            yMax = min(yMin+rangeY-1,imageSizeY);
            yMin = max(yMax-rangeY+1,1);
            zMin = max(iZm-((rangeZ-1)/2),1);
            zMax = min(zMin+rangeZ-1,imageSizeZ);
            zMin = max(zMax-rangeZ+1,1);
            
            % get the local image from the cube defined by the limits
            imageLocal = imageLast5(xMin:xMax,yMin:yMax,zMin:zMax,:);
            
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
            
            % get range of bgMean and bgStd to update
            
            
            
            %put values in matrix representing image
            %             bgMean(iX,iY,iZ) = bgMean1;
            %             bgStd(iX,iY,iZ) = bgStd1;
            bgMean(xMin:xMax,yMin:yMax,zMin:zMax) = bgMean1;
            bgStd(xMin:xMax,yMin:yMax,zMin:zMax) = bgStd1;
            
            timeCount = timeCount + 1;
            %             progressText(timeCount/timeTotal,['Calculating local background at {x,y,z}={',int2str(iX),',',int2str(iY),',',int2str(iZ),'}']);
            if verbose
                progressText(timeCount/timeTotal,sprintf('Calculating local background at {x,y,z}={%1.3d,%1.3d,%1.3d}',iX,iY,iZ));
            end
        end
    end
end
