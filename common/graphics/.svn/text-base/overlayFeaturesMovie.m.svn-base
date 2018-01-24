function overlayFeaturesMovie(movieInfo,startend,saveMovie,movieName,...
    filterSigma,showRaw)
%Overlays detected features obtained via detectSubResFeatures2D_Movie on movies
%
%SYNPOSIS overlayFeaturesMovie(movieInfo,startend,saveMovie,movieName,filterSigma)
%
%INPUT  movieInfo   : Output of detectSubResFeatures2D_Movie.
%       startend    : Row vector indicating first and last frame to
%                     include in movie. Format: [startframe endframe]. 
%                     Optional. Default: [1 (maximum available frame)]            
%       saveMovie   : 1 to save movie (as Quicktime), 0 otherwise.
%                     Optional. Default: 0
%       movieName   : filename for saving movie. 
%                     Optional. Default: FeaturesMovie (if saveMovie = 1).
%       filterSigma : 0 to overlay on raw image, PSF sigma to overlay on image
%                     filtered with given filterSigma. 
%                     Optional. Default: 0
%       showRaw     : 1 to add raw movie to the left of the movie with
%                     tracks overlaid, 2 to add raw movie at the top of
%                     the movie with tracks overlaid, 0 otherwise.
%                     Optional. Default: 0.
%
%OUTPUT If movie is to be saved, the QT movie is written into directory
%       where TIFFs are located
%
%Khuloud Jaqaman, August 2007

%% input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--overlayFeaturesMovie: Incorrect number of input arguments!');
    return
end

%record directory before start of function
startDir = pwd;

%ask user for images
[fName,dirName] = uigetfile('*.tif','specify first image in the stack - specify very first image, even if not to be plotted');

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileList=getFileStackNames([dirName,fName]);
    numFrames = length(outFileList);
    
    %read first image to get image size
    currentImage = imread(outFileList{1});
    [isx,isy] = size(currentImage);

else %else, exit
    
    disp('--overlayTracksMovieNew: Bad file selection');
    return

end

%check startend and assign default if necessary
if nargin < 2 || isempty(startend)
    startend = [1 numFrames];
else
    startend(2) = min(startend(2),numFrames); %make sure that last frame does not exceed real last frame
end

%check whether to save movie
if nargin < 3 || isempty(saveMovie)
    saveMovie = 0;
end

%check name for saving movie
if saveMovie && (nargin < 4 || isempty(movieName))
    movieName = 'featuresMovie.mov';
end

%check whether to use filtered images
if nargin < 5 || isempty(filterSigma)
    filterSigma = 0;
end

%check whether to put raw movie adjacent to movie with tracks overlaid
if nargin < 6 || isempty(showRaw)
    showRaw = 0;
end

%keep only the frames of interest
outFileList = outFileList(startend(1):startend(2));

%read specified images into ImageStack
ImageStack = zeros(isx,isy,length(outFileList));
h = waitbar(0,'reading images');
for i=1:length(outFileList)
    currentImage = imread(outFileList{i});
    ImageStack(:,:,i) = currentImage;
    if (mod(i,10)==0), waitbar(i/length(outFileList)); end
end
close(h);

%filter images if requested
if filterSigma
    for i = 1 : length(outFileList)
        imageStack(:,:,i) = Gauss2D(imageStack(:,:,i),filterSigma);
    end
end

%initialize QT movie if it is to be saved
if saveMovie
    evalString = ['MakeQTMovie start ',movieName];
    eval(evalString);
end

%retain only the movieInfo of the frames of interest
if isempty(movieInfo)
    movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),...
        startend(2)-startend(1)+1,1);
else
    movieInfo = movieInfo(startend(1):startend(2));
end

%get image size
[imSizeX,imSizeY] = size(ImageStack(:,:,1));
imageRange = [1 imSizeX; 1 imSizeY];


%% make movie

%go to directory where movie will be saved
cd(dirName);

%go over all specified frames
for iFrame = 1 : length(movieInfo)
    
    %plot image in current frame
    clf;
    
    switch showRaw
        case 1
            axes('Position',[0 0 0.495 1]);
            imshow(ImageStack(:,:,iFrame),[]);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
                textDeltaCoord,num2str(iFrame),'Color','white');
            axes('Position',[0.505 0 0.495 1]);
            imshow(ImageStack(:,:,iFrame),[]);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
        case 2
            axes('Position',[0 0.505 1 0.495]);
            imshow(ImageStack(:,:,iFrame),[]);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
                textDeltaCoord,num2str(iFrame),'Color','white');
            axes('Position',[0 0 1 0.495]);
            imshow(ImageStack(:,:,iFrame),[]);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
        otherwise
            axes('Position',[0 0 1 1]);
            imshow(ImageStack(:,:,iFrame),[]);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
                textDeltaCoord,num2str(iFrame),'Color','white');
    end
    
    %plot features
    if ~isempty(movieInfo(iFrame).xCoord)
        plot(movieInfo(iFrame).xCoord(:,1),movieInfo(iFrame).yCoord(:,1),'ro','MarkerSize',4);
    end
        
    %add frame to movie if movie is saved
    if saveMovie
        MakeQTMovie addfigure
    end

    %pause for a moment to see frame
    pause(0.1);

end

%finish movie
if saveMovie==1
    MakeQTMovie finish
end

%% change directory back to original
cd(startDir);

%% ~~~ end ~~~