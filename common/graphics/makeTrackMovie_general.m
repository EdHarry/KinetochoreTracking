function []=makeTrackMovie_general(trackInfo,axisVector,startend,dragtailLength,recordingVariable,moviename)
% makes a tracking movie with dragtails using the original tiff images and 
% the tracking info matrix (as produced by Khuloud's tracker)
% SYNPOSIS: []=makeTrackMovie_general(trackInfo,axisVector,startendVector,dragtailLength,recordingVariable,moviename)
%
% INPUT:    trackInfo      = trackInfo matrix from Khuloud's tracker
%           axisVector(opt)= vector containing the axis specifications for
%               the image in the format [xstart xend ystart yend]; 
%   `           default value is: [1 imageSizeX 1 imageSizeY]
%           startend (opt) = vector containing desired first and last frame
%               of the movie in the format [startframe endframe]; 
%               default value is: [1 (maximum available frame)]            
%           dragtailLength (opt) = length of drag tail; default: 10 frames
%               NOTE: if dragtailLength=0, then no dragtail is made; if the 
%               dragtail should have maximum length (as long as the movie),
%               then set it to any value longer than the movie, e.g. 10000 
%           recordingVariable (opt)   = should the movie be written into a 
%               Quicktime movie, value 1/0 corresponds to yes/no
%               default value is: 0 (no recording)
%           moviename(opt) = filename for movie, if recordingVariable=1
%
% OUTPUT: none, but if movie is recorded, then QT movie is written into
%           directory where TIFFs are located
%
%         In the Quicktime movie, the overlaid markers mean the following:
%           round markers:    detected objects
%           stars:            closed gaps
%         The colors mean the following:
%           *round markers*:
%           magenta = objects appear OR disappear (they are present in
%                     either the first or the last frame of the movie)
%           white = objects are stationary (they are present for the
%                   entire course of the movie, including the first and
%                   last frame)
%           red = objects appear AND disappear (they are not present in the
%                 first or last frame), which makes it possible to
%                 determine their exact lifetime
%           green = object appearance (first frame of trajectory)
%           yellow = object disappearance (last frame of trajectory)
%           *stars*:
%           blue = 'bad gap' (gap is longer than adjacent trajectory pieces)
%           cyan = 'good gap' (gap is shorter than adjacent trajectory pieces)
%
% NOTE1: this version of the function does not require the presence of the
% images loaded into the workspace, rather, the images are loaded
% temporarily from the location specified where they are stored
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: dloerke
% DATE: 25-Jan-2007
%
% last modified 01/26/2006
%

%% =======================================================================
%%======        Step 0: initialize and set default values
%% =======================================================================

% size of trackinfo matrix
[tix,tiy] = size(trackInfo);
tinumf = round(tiy/8);


% set optional parameters to default values if necessary

% start and end frames
if nargin>2 && ~isempty(startend)
    startendVector = startend;   
else
    startendVector = [1 tinumf];
end

% dragtail length
if nargin>3 && ~isempty(dragtailLength)
    dtLength = dragtailLength;
else
    dtLength = 10;
end

% recording status
if nargin>4 && ~isempty(recordingVariable)
    recStatus = recordingVariable;
    % if recording is desired, initialize the quicktime movie
    if (recStatus==1) 
        % if movie name is specified, use name, else call it 'Trackmovie'
        if nargin>5 && ~isempty(moviename)
            evalString = ['MakeQTMovie start ',moviename];
        else
            evalString = ('MakeQTMovie start Trackmovie');
        end
        eval(evalString);
    end
else
    recStatus = 0;
end



%% =======================================================================
%%======        Step 1: read images in specified directory
%% =======================================================================

% record directory before start of function
startDir = pwd;


% determine image file names 
[fName, dirName] = uigetfile('*.tif','specify first image in the stack');
% NOTE: here, the user should specify the first available image file in the 
% stack, such as 'imageName_001', regardless of whether it will be used. If
% the desired first frame in the movie (specified in the input vector
% startend) is different from 1, then the first image(s) will be discarded
% in the subsequent lines


if(isa(fName,'char') & isa(dirName,'char'))
     
    % Recover all file names from the stack
    outFileList=getFileStackNames([dirName,fName]);
    
    % crop the portion of interest from the filename list, BUT: don't
    % read more images than there are frames in the trackInfo matrix
    startendVector(2) = min(length(outFileList),startendVector(2));
    outFileList=outFileList(startendVector(1):startendVector(2));
        
    % read specified images into ImageStack
    h= waitbar(0,'reading images');
    for i=1:length(outFileList)
        currentImage = imread(outFileList{i});
        ImageStack(:,:,i) = currentImage;
        if (mod(i,10)==0), waitbar(i/length(outFileList)); end
    end
    close(h);
end

% if axisVector isn't specified - either if there's no entry, or if the 
% input is [] - then set axisVector to the image size
[isx,isy]=size(currentImage);
if nargin<2
    axisVector = [1 isy 1 isx];
elseif length(axisVector)==0
    axisVector = [1 isy 1 isx];
end



%% =======================================================================
%% = Step 2: from trackInfo matrix, determine lifetimes and status matrices
%% =======================================================================

startp = startendVector(1); 
endp = startendVector(2); 
startindex = 8*(startp-1)+1;
endindex = 8*endp;

[Mlft,Mstat,Mx,My,Mda]=findLifetimesStatusSimple(trackInfo(:,startindex:endindex));



%% =======================================================================
%% ======       Step 3: make movie with the results
%% =======================================================================

cd(dirName);

makeMovieFromLFTstatus_general(ImageStack,Mlft,Mstat,Mx,My,Mda,dtLength,axisVector,recStatus);

% if movie is recorded, finish movie after this step
if recStatus==1
    MakeQTMovie finish
end


%% change directory back to original
cd(startDir);


end % of function



%% =======================================================================
%% =======================================================================
%% ====
%% ====                          SUBFUNCTIONS
%% ==== 
%% =======================================================================
%% =======================================================================


function []=makeMovieFromLFTstatus_general(imageStack,lftTemp,statTemp,xTemp,yTemp,daTemp,dtLength,axisVector,recinput);

% check if writing to a QT movie is desired
if nargin>8
    if recinput==1
        rec=1;
    else
        rec=0;
    end
else
    rec=0;
end

% number of available frames is determined from size of image stack
[isx,isy,nf]=size(imageStack);


for i=1:nf
    % vectors with lifetimes, status and position in current frame
    tvec_lft = full(lftTemp(:,i));
    tvec_stat = full(statTemp(:,i));
    tvec_da = full(daTemp(:,i));
    tvec_x = full(xTemp(:,i));
    tvec_y = full(yTemp(:,i));
    
    % positions of trajectories of different categories
    pstat1 = find(tvec_stat==1);
    pstat2 = find(tvec_stat==2);
    pstat3 = find(tvec_stat==3);
    % positions of gaps of different categories
    pstat4 = find(tvec_stat==4);
    pstat5 = find(tvec_stat==5);
    % positions of appearance and disappearances
    astat = find(tvec_da==1);
    dastat = find(tvec_da==-1);

    
    % plot image:
    % for better contrast, the [low high] settings are not set to default, 
    % but to the min and max values of the portion of the image that is 
    % actually displayed
    useImage = imageStack(:,:,i);
    cropImage = useImage(axisVector(3):axisVector(4),axisVector(1):axisVector(2));
    lo = min(cropImage(:));
    hi = max(cropImage(:));
    
    clf;
    imshow(imageStack(:,:,i),[lo hi]); hold on; 
    axis(axisVector);
    
    % make dragtails
    % NOTE: in the present form, dragtails are only shown for existing 
    % objects; once an object disappears, its past trajectory/dragtail 
    % will disappear as well, but this can of course be changed according
    % to the needs of the user
    if i>1
        % positions of currently existing objects 
        xcurr = full(xTemp(:,i));
        ycurr = full(yTemp(:,i));
        xycurrpos = find((xcurr>0)&(ycurr>0));
        xcurrUse = xcurr(xycurrpos);
        ycurrUse = ycurr(xycurrpos);
        
        xComplete = xcurrUse;
        yComplete = ycurrUse;
        
        % if a dragtail is specified (dtLength>0), then determine the 
        % previous positions (of the currently existing objects)
        % along the specified length of the dragtail and plot them
        if dtLength>0
            pmin = max(1,i-dtLength);
            for p=(i-1):-1:pmin
                xprev = full(xTemp(:,p));
                yprev = full(yTemp(:,p));
                xprevUse = xprev(xycurrpos);
                yprevUse = yprev(xycurrpos);
                xprevUse(xprevUse==0)=nan;
                yprevUse(yprevUse==0)=nan;
                xComplete = [xprevUse xComplete];
                yComplete = [yprevUse yComplete];
            end % of for
            
            % now plot all individual dragtails
            for k=1:length(xycurrpos)
                plot(xComplete(k,:),yComplete(k,:),'r-');
            end
        end % of if dragtail is >0
        
        
    end
    
    % plot the positions of the current objects, where:
    % round markers:    detected objects
    % stars:            closed gaps
    
    % magenta = status 2, objects appear OR disappear
    plot(tvec_x(pstat2),tvec_y(pstat2),'mo','MarkerSize',4,'LineWidth',1);
    % white = status 3, objects are stationary (present for the entire movie)
    plot(tvec_x(pstat3),tvec_y(pstat3),'wo','MarkerSize',4,'LineWidth',1);
    % red = status 1, objects appear AND disappear
    plot(tvec_x(pstat1),tvec_y(pstat1),'ro','MarkerSize',4,'LineWidth',1);
    
    % green = astat, object appears in this frame
    plot(tvec_x(astat),tvec_y(astat),'go','MarkerSize',4,'LineWidth',1);
    % yellow = dastat, object disappears in this frame
    plot(tvec_x(dastat),tvec_y(dastat),'yo','MarkerSize',4,'LineWidth',1);
    
    % blue = 'bad gap' (gap is longer than adjacent trajectory pieces)
    plot(tvec_x(pstat5),tvec_y(pstat5),'b*','MarkerSize',6,'LineWidth',1);
    % cyan = 'good gap' (gap is shorter than adjacent trajectory pieces)
    plot(tvec_x(pstat4),tvec_y(pstat4),'c*','MarkerSize',6,'LineWidth',1);
    
    text(axisVector(1)+10,axisVector(3)+10,num2str(i),'Color','white');
        
    % figure is recorded if indicated
    if rec==1
        MakeQTMovie addfigure
    end
    % pause to enjoy the movie    
    pause(0.1);
end % of for

end % of sub function




