function [averageImage,sdImage,hotImage,xPattern, yPattern, xPatternL, yPatternL]=cameraEvaluation(darkFrames, testDynamicPattern, nImagesForCorr, nLags, imgName)
%CAMERAEVALUATION tests cameras for static and dynamic patterns
%
% The code performs an analysis on a series of dark images. 
%       First, the average and standard deviation of the image series are
%       calculated. These are the static patterns. For better visibility,
%       hot pixels in the static patterns are masked with the robust mean
%       in the images.
%       To find dynamic patterns, the dark images are corrected by the
%       average and divided by the standard deviation. Then, every single
%       image is analyzed twice for dynamic correlation: Once by looking
%       for correlation on the entire corrected images, and once on images
%       with a 20% border removed. For every lag, the number of significant
%       correlations is counted (95% level). If there are more than 5%
%       significant lags in any given direction (=above the red line in the
%       graphs), there is a correlation.
%       The code writes the number of hot pixels to the commandline.
%
% SYNOPSIS [averageImage,sdImage,hotImage,xPattern, ...
%               yPattern, xPattern2, yPattern2]= ...
%               cameraEvaluation(darkFrames, ...
%               testDynamicPattern, nImagesForCorr, nLags) 
%
% INPUT    darkFrames : 3D, 4D or 5D image stack (will be converted to a 3D
%                       image stack) 
%          testDynamicPattern : (opt) whether to test for a dynamic pattern
%                       (significant autocorrelation left after the average
%                       image has been subtracted and the standard
%                       deviation has been normalized). [0/{1}]
%          nImagesForCorr : how many images should be chosen to run the
%                       autocorrelation analysis. inf = all. Default:
%                       min(25,nImages). Careful: Increasing this number
%                       will slow down things tremendously.
%          nLags :      Determines the maximum lag for the autoCorrelation
%                       as nLags*max(imageSize). Default: 3.1
%          imgName :    Name of image
%
% OUTPUT   averageImage : average of the image series
%          sdImage      : standard deviation of the image series
%          hotImage     : binary image of hot pixels
%          xPattern, yPattern, xPattern2, yPattern2
%                       : Logical arrays with the incidence of significant
%                         correlations.
%
% c: jonas, 10/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TEST INPUT

% default parameters
%def_borderPercent = 20;
def_testDynamicPattern = 1;
def_nImagesForCorr = 25;
def_nLags = 3.1;

if nargin < 2 || isempty(testDynamicPattern)
    testDynamicPattern = def_testDynamicPattern;
end
if nargin < 3 || isempty(nImagesForCorr)
    nImagesForCorr = def_nImagesForCorr;
end
if nargin < 4 || isempty(nLags)
    nLags = def_nLags;
end
if nargin < 5 || isempty(imgName)
    imgName = inputname(1);
end

% convert image
movieSize = size(darkFrames);
darkFrames = reshape(darkFrames,movieSize(1),movieSize(2),[]);
movieSize = size(darkFrames);

% take average and std
averageImage = mean(darkFrames,3);
sdImage = std(darkFrames,0,3);

[robMean] = robustMean(averageImage(:));
[robStd, dummy, inlierIdx] = robustMean(sdImage(:));

hotImage = averageImage > (robMean + 10 * robStd/sqrt(movieSize(3)));
numberOfHotPixels = nnz(hotImage);
ratioOfHotPixels = numberOfHotPixels/prod(movieSize(1:2));


% correct averageImage, stdImage for hot pixels
averageImageC = averageImage;
averageImageC(hotImage) = robMean;
sdImageC = robStd * ones(movieSize(1:2));
sdImageC(inlierIdx) = sdImage(inlierIdx);


figure('Name',sprintf('Average background %s',imgName))
imshow(averageImageC,[]);
set(gcf,'Name','Average background');
colormap('jet')


figure('Name',sprintf('Std background %s',imgName))
imshow(sdImageC,[]);
colormap('jet')

figure('Name',sprintf('Average background (no hot pixels) %s',imgName)),
histogram(averageImageC(~hotImage));
figure('Name',sprintf('Std background (no hot pixels) %s',imgName));
histogram(sdImageC(~hotImage));

% subtract static pattern, divide by sd
darkFrames = darkFrames - repmat(averageImage,[1,1,movieSize(3)]);
darkFrames = darkFrames ./ repmat(sdImage,[1,1,movieSize(3)]);

% plot average corrected darkFrame
% correctedAverageImage = mean(darkFrames,3);
% correctedSdImage = std(darkFrames,0,3);
% uiViewPanel,imshow(correctedAverageImage,[])
% set(gcf,'Name','average of corrected background')
% colormap('jet')
% uiViewPanel,imshow(correctedSdImage,[])
% set(gcf,'Name','std of corrected background')
% colormap('jet')
% 
% figure('Name','Average corrected background')
% histogram(correctedAverageImage(:));
% 
% figure('Name','Std corrected background')
% histogram(correctedSdImage(:));

if ~testDynamicPattern
    
    % assign empty output
    [xPattern, yPattern, xPatternL, yPatternL] = deal([]);
    
else % test the dynamic pattern

% get dynamic pattern. Do either 10 times the number of pixels along the
% longer side of the image, or half the total number of pixels in the image
%nCorr = min(max(10*movieSize(1:2)),ceil(prod(movieSize(1:2))/2));

% correlate 3 times longer side
% run for 25 random images
imgList = randperm(movieSize(3));

imgList = imgList(1:nImagesForCorr);
nCorr = ceil(nLags*max(movieSize(1:2)));
[xPattern, yPattern] = deal(zeros(nCorr+1,nImagesForCorr));
progressText(0,'Correlations') % Create text
ct = 1;
for i = imgList
    currentFrame = darkFrames(:,:,i);
    out = autoCorr(currentFrame(:),nCorr,-1);
    xPattern(:,ct) = out(:,1);
    currentFrame = currentFrame';
    out = autoCorr(currentFrame(:),nCorr,-1);
    yPattern(:,ct) = out(:,1);
    progressText(ct/nImagesForCorr);
    ct = ct+1;
end

% check for significance: is autocorrelation larger than 
% 1.96/sqrt(numberOfPixels)?
numberOfPixels = prod(movieSize(1:2));
% significance threshold: choose such that only 1 pixel should be
% significant
thresh = norminv(1-0.5/numberOfPixels)/sqrt(numberOfPixels);
xPatternL = abs(xPattern) > thresh;
yPatternL = abs(yPattern) > thresh;

% plot

% correlations
figure('Name',sprintf('average correlation rows (npts=%i, nimg=%i) %s',numberOfPixels,nImagesForCorr, imgName))
m=mean(xPattern,2);
s=std(xPattern,0,2);
area(0:nCorr,m+s,'BaseValue',-1,'LineStyle','none','FaceColor',[0.8,0.8,0.8])
hold on
area(0:nCorr,m-s,'BaseValue',-1,'LineStyle','none','FaceColor','w')
plot(0:nCorr,m)
plot(0:nCorr,repmat(thresh,nCorr+1,1),'r')
xlabel(sprintf('Pixel lags (%i pix/row)',movieSize(1)))
ylabel('Normalized autocorrelation')

figure('Name',sprintf('average correlation rows (npts=%i, nimg=%i) %s',numberOfPixels,nImagesForCorr,imgName))
m=mean(yPattern,2);
s=std(yPattern,0,2);
area(0:nCorr,m+s,'BaseValue',-1,'LineStyle','none','FaceColor',[0.8,0.8,0.8])
hold on
area(0:nCorr,m-s,'BaseValue',-1,'LineStyle','none','FaceColor','w')
plot(0:nCorr,m)
plot(0:nCorr,repmat(thresh,nCorr+1,1),'r')
xlabel(sprintf('Pixel lags (%i pix/row)',movieSize(2)))
ylabel('Normalized autocorrelation')

% sum of significant correlations
figure('Name',sprintf('Number of significant correlations in rows (npts=%i, nimg=%i) %s',numberOfPixels,nImagesForCorr,imgName))
plot(0:nCorr,sum(xPatternL,2),'.')
xlabel(sprintf('Pixel lags (%i pix/row)',movieSize(2)))
ylabel(sprintf('Number of significant correlations out of %i',nImagesForCorr))
ylim([-0.2,nImagesForCorr+0.2])
grid on

figure('Name',sprintf('Number of significant correlations in rows (npts=%i, nimg=%i) %s',numberOfPixels,nImagesForCorr,imgName))
plot(0:nCorr,sum(yPatternL,2),'.')
xlabel(sprintf('Pixel lags (%i pix/row)',movieSize(2)))
ylabel(sprintf('Number of significant correlations (max %i)',nImagesForCorr))
ylim([-0.2,nImagesForCorr+0.2])
grid on

% % the first row is obviously not informative
% xPatternL(1,:) = false;
% yPatternL(1,:) = false;
% 
% % uiViewPanel,imshow(xPattern',[])
% % uiViewPanel,imshow(yPattern',[])
% 
% % sum the incidence of significant lags. Normalize by the number of frames
% sx = sum(xPatternL,2)/movieSize(3);
% sy = sum(yPatternL,2)/movieSize(3);
% 
% figure('Name','Correlation X'),stem([0:nCorr-1]',sx);
% hold on
% plot([0:nCorr-1],0.05*ones(nCorr,1),'r')
% figure('Name','Correlation Y'),stem([0:nCorr-1],sy);
% hold on
% plot([0:nCorr-1],0.05*ones(nCorr,1),'r')

% % do again, but without 20% pixels all around
% xPercent = floor(movieSize(1:2)/(borderPercent/100));
% if all(xPercent > 2)
%     
%     darkFrames = darkFrames(xPercent(1)+1:end-xPercent(1),...
%         xPercent(2)+1:end-xPercent(2),:);
%     movieSize = size(darkFrames);
%     nCorr = min(max(10*movieSize(1:2)),ceil(prod(movieSize(1:2))/2));
%     [xPattern2, yPattern2] = deal(zeros(nCorr,movieSize(3)));
%     tic
%     for i = 1:movieSize(3)
%         [xPattern2(:,i), yPattern2(:,i)] = ...
%             normACfunc(darkFrames(:,:,i), -1);
%         disp(sprintf('iteration: %i/%i, time: %9.3f',i,movieSize(3),toc))
%     end
% 
%     % check for significance: is autocorrelation larger than
%     % 1.96/sqrt(numberOfPixels)?
%     numberOfPixels = prod(movieSize(1:2));
%     xPattern2 = abs(xPattern2) > 1.96/sqrt(numberOfPixels);
%     yPattern2 = abs(yPattern2) > 1.96/sqrt(numberOfPixels);
% 
%     % the first row is obviously not informative
%     xPattern2(1,:) = logical(0);
%     yPattern2(1,:) = logical(0);
% 
%     % uiViewPanel,imshow(xPattern',[])
%     % uiViewPanel,imshow(yPattern',[])
% 
%     % sum the incidence of significant lags. Normalize by the number of frames
%     sx = sum(xPattern2,2)/movieSize(3);
%     sy = sum(yPattern2,2)/movieSize(3);
% 
%     figure('Name',...
%         sprintf('Correlation X (minus 2x %i border pix)',xPercent(1)));
%     stem([0:nCorr-1]',sx);
%     hold on
%     plot([0:nCorr-1],0.05*ones(nCorr,1),'r')
%     figure('Name',...
%         sprintf('Correlation Y (minus 2x %i border pix)',xPercent(2)));
%     stem([0:nCorr-1],sy);
%     hold on
%     plot([0:nCorr-1],0.05*ones(nCorr,1),'r')
% end
% 
 end % if testDynamicPattern

% display number of hot pixels
disp(sprintf('Number of hot pixels: %i (%2.3f%%)',numberOfHotPixels, ...
    100*ratioOfHotPixels));

