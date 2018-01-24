function [intensityHist,numPeaks,firstGaussMean,firstGaussStd,ratioPeak2toPeak1,...
    ratioPeak3toPeak2,numFeatures,errFlag] = singleFluoroAnalysis(movieInfo,...
    movieName,startFrame,endFrame,alpha,variableMean,variableStd,maxNumGauss,...
    plotResults)
%SINGLEFLUOROANALYSIS looks for peaks in the intensity histogram in each frame and plots their characteristics.
%
%SYNOPSIS [intensityHist,numPeaks,firstGaussMean,firstGaussStd,ratioPeak2toPeak1,...
%    ratioPeak3toPeak2,numFeatures,errFlag] = singleFluoroAnalysis(movieInfo,...
%    movieName,startFrame,endFrame,alpha,variableMean,variableStd,maxNumGauss,...
%    plotResults)
%
%INPUT  
%   Mandatory
%       movieInfo    : Array of size equal to the number of frames
%                      in movie, containing the fields:
%             .xCoord      : Image coordinate system x-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%             .yCoord      : Image coordinate system y-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%                            Optional. Skipped if problem is 1D. Default: zeros.
%             .zCoord      : Image coordinate system z-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%                            Optional. Skipped if problem is 1D or 2D. Default: zeros.
%             .amp         : Amplitudes of PSFs fitting detected features. 
%                            1st column for values and 2nd column 
%                            for standard deviations.
%   Optional
%       movieName    : Name of movie. Needed to put as title of plot if plotting.
%                      Default: 'movie1'.
%       startFrame   : Frame number at which to start analysis.
%                      Default: 1.
%       endFrame     : Frame number at which to end analysis.
%                      Default: Last frame in movie.
%       alpha        : Alpha-value for the statistical test that compares the
%                      fit of n+1 Gaussians to the fit of n Gaussians.
%                      Default: 0.05.
%       variableMean : 0 if assuming the fixed relationship
%                      (mean of nth Gaussian) = n * (mean of 1st Gaussian).
%                      1 if there is no relationship between the means of
%                      different Gaussians.
%                      Default: 0.
%       variableStd  : 0 if assuming that all Gaussians have the same
%                      standard deviation. 1 if there is no relationship 
%                      between the standard deviations of different
%                      Gaussians, 2 if assuming the relationship 
%                      (std of nth Gaussian) = sqrt(n) * (std of 1st Gaussian). 
%                      variableStd can equal 2 only if variableMean is 0.
%                      Default: 0.
%       maxNumGauss  : upper limit to the number of fitted Gaussians.
%                      Default: 10.
%       plotResults  : 1 if results are to be plotted, 0 otherwise.
%                      Default: 1.
%
%OUTPUT intensityHist    : Array of structures with field gaussParam as outputed by
%                          fitHistWithGaussians for each analyzed frame.
%       numPeaks         : Number of intensity peaks detected in each
%                          analyzed frame.
%       firstGaussMean   : Mean of first fitted Gaussian in each analyzed
%                          frame (i.e. position of first peak in intensity
%                          histogram).
%       firstGaussStd    : Standard deviation of first fitted Gaussian in
%                          each analyzed frame (which is proportional to 
%                          the width of the first peak in intensity histogram).
%       ratioPeak2toPeak1: Ratio of amplitude of second peak to first peak
%                          per analyzed frame. NaN indicates that a frame 
%                          does not have a second peak.
%       ratioPeak3toPeak2: Ratio of amplitude of third peak to second peak
%                          per analyzed frame. NaN indicates that a frame 
%                          does not have a second peak.
%       numFeatures      : Number of detected features in each analyzed
%                          frame. NaN indicates frame hasn't been analyszed.
%       errFlag          : 0 if function executes normally, 1 otherwise;
%
%Khuloud Jaqaman, April 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--analyzeIntensities: Function needs at least 1 input argument!');
    return
end

%find number of frames in movie
numFrames = length(movieInfo);

%assign default values of optional input variables
movieName_def = 'movie1';
startFrame_def = 1;
endFrame_def = numFrames;
alpha_def = 0.05;
variableMean_def = 0;
variableStd_def = 0;
maxNumGauss_def = 10;
plotResults_def = 1;

%check movieName
if nargin < 2 || isempty(movieName)
    movieName = movieName_def;
end    

%check startFrame
if nargin < 3 || isempty(startFrame)
    startFrame = startFrame_def;
else
    if startFrame < 1 || startFrame > numFrames
        disp('--analyzeIntensities: "startFrame" should be between 1 and number of frames in movie!');
        errFlag = 1;
    end
end

%check endFrame
if nargin < 4 || isempty(endFrame)
    endFrame = endFrame_def;
else
    if endFrame < startFrame || endFrame > numFrames
        disp('--analyzeIntensities: "endFrame" should be between "startFrame" and number of frames in movie!');
        errFlag = 1;
    end
end

%check alpha
if nargin < 5 || isempty(alpha)
    alpha = alpha_def;
else
    if alpha < 0 || alpha > 1
        disp('--analyzeIntensities: "alpha" should be between 0 and 1!');
        errFlag = 1;
    end

end

%check variableMean
if nargin < 6 || isempty(variableMean)
    variableMean = variableMean_def;
else
    if ~any(variableMean == [0,1])
        disp('--analyzeIntensities: "variableMean" should be 0 or 1!');
        errFlag = 1;
    end
end

%check variableStd
if nargin < 7 || isempty(variableStd)
    variableStd = variableStd_def;
else
    if ~any(variableStd == [0,1,2])
        disp('--analyzeIntensities: "variableStd" should be 0, 1 or 2!');
        errFlag = 1;
    end
end

%check maxNumGauss
if nargin < 8 || isempty(maxNumGauss)
    maxNumGauss = maxNumGauss_def;
else
    if maxNumGauss < 1
        disp('--analyzeIntensities: "maxNumGauss" should be at least 1!');
        errFlag = 1;
    end
end

%check plotResults
if nargin < 9 || isempty(plotResults)
    plotResults = plotResults_def;
else
    if ~any(plotResults == [0,1])
        disp('--analyzeIntensities: "plotResults" should be 0 or 1!');
        errFlag = 1;
    end
end

%exit if there are problem in input variables
if errFlag
    disp('--analyzeIntensities: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Intensity analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of features in each frame
numFeatures = NaN*ones(1,numFrames);
for i = startFrame : endFrame
    if ~isempty(movieInfo(i).amp)
        numFeatures(i) = length(movieInfo(i).amp(:,1));
    else
        numFeatures(i) = 0;
    end
end

%fit the amplitude distribution in each frame with Gaussians
%don't consider frame if it has less than 10 features
intensityHist(1:numFrames) = struct('gaussParam',[]);
for i = startFrame : endFrame
    if numFeatures(i) >= 10
        [numObsPerBin,binCenter,gaussParam,errFlag] = ...
            fitHistWithGaussians(movieInfo(i).amp(:,1),alpha,variableMean,...
            variableStd,0,maxNumGauss,2);
        intensityHist(i).gaussParam = gaussParam;
    end
end

%get number of peaks found in each frame
numPeaks = NaN*ones(1,numFrames);
for i = startFrame : endFrame
    if ~isempty(intensityHist(i).gaussParam)
        numPeaks(i) = size(intensityHist(i).gaussParam,1);
    end
end

%get mean and std of first Gaussian in each frame
firstGaussMean = NaN*ones(1,numFrames);
firstGaussStd = NaN*ones(1,numFrames);
for i = startFrame : endFrame
    if ~isempty(intensityHist(i).gaussParam)
        firstGaussMean(i) = intensityHist(i).gaussParam(1,1);
        firstGaussStd(i) = intensityHist(i).gaussParam(1,2);
    end
end

%get ratio of amplitude of 2nd peak to 1st peak, and 3rd peak to 2nd peak
ratioPeak2toPeak1 = NaN*ones(1,numFrames);
ratioPeak3toPeak2 = NaN*ones(1,numFrames);
for i = startFrame : endFrame
    if ~isempty(intensityHist(i).gaussParam)
        switch numPeaks(i)
            case 2 %if there are two peaks
                ratioPeak2toPeak1(i) = ...
                    (intensityHist(i).gaussParam(2,3)/intensityHist(i).gaussParam(2,2))/...
                    (intensityHist(i).gaussParam(1,3)/intensityHist(i).gaussParam(1,2));
            case 3 %if there are three peaks
                ratioPeak2toPeak1(i) = ...
                    (intensityHist(i).gaussParam(2,3)/intensityHist(i).gaussParam(2,2))/...
                    (intensityHist(i).gaussParam(1,3)/intensityHist(i).gaussParam(1,2));
                ratioPeak3toPeak2(i) = ...
                    (intensityHist(i).gaussParam(3,3)/intensityHist(i).gaussParam(3,2))/...
                    (intensityHist(i).gaussParam(2,3)/intensityHist(i).gaussParam(2,2));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot results
if plotResults
    
    %open new figure window
    figure('Name',movieName,'NumberTitle','off');

    %plot number of features per frame in 1st sub-plot
    subplot(4,1,1)
    hold on;
    if any(~isnan(numFeatures))
        plot(numFeatures,'k');
        axis([0 numFrames 0 1.1*max(numFeatures)]);
    end
    title(['number of features for ' movieName]);

    %plot number of peaks per frame in 2nd sub-plot
    subplot(4,1,2)
    hold on;
    if any(~isnan(numPeaks))
        plot(numPeaks,'k')
        axis([0 numFrames 0 max(numPeaks)+1]);
    end
    title('number of peaks');
    
    %plot mean and std of first fitted Gaussian in 3rd sub-plot
    subplot(4,1,3)
    hold on
    if any(~isnan(firstGaussMean))
        plot(firstGaussMean,'k')
        plot(firstGaussStd,'r')
        axis([0 numFrames 0 1.1*max([firstGaussMean firstGaussStd])]);
    end
    title('mean (black) and std (red) of 1st peak');
    
    %plot ratio of amplitudes of 2nd peak to 1st peak and 3rd peak to 2nd
    %peak
    subplot(4,1,4)
    hold on
    if any(~isnan(ratioPeak2toPeak1))
        plot(ratioPeak2toPeak1,'k','marker','.')
        plot(ratioPeak3toPeak2,'r','marker','.')
        axis([0 numFrames 0 1.1*max([ratioPeak2toPeak1 ratioPeak3toPeak2])]);
    end
    title('ratio of amplitudes of 2nd to 1st peak (black) and 3rd to 2nd peak (red)');
    xlabel('frame number');
    
end %(if plotRes)


%%%%% ~~ the end ~~ %%%%%
