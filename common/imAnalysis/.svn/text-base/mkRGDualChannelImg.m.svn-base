function mkRGDualChannelImg(varargin)
%mkRGDualChannelImg: This function creates Red/Green duo-channel images.
%
% SYNOPSIS:
%    mkRGDualChannelImg;
%    mkRGDualChannelImg('outFileName',outFileName,'redImgIScale',redImgIScale,'greenImgIScale',greenImgIScale);
%
% OPTIONAL PAR/VALUE INPUT:
%        PAR
%    -------------
%    'outFileName': The common prefix name for the output files. Default is the
%                   same name as the image file name of the first channel with 'rg_' as prefix.
%    'redImgIScale': (red channel) A vector of two numbers between 0 and 1 that specifies
%                    the percentage of cut-offs at the lower and upper end of image intensity.
%    'greenImgIScale': (green channel) A vector of two numbers between 0 and 1 that specifies
%                    the percentage of cut-offs at the lower and upper end of image intensity.
%    'greenToRedIRatio': An intensity ratio of green over red can also be
%                    specified so that signals in one channel may be enhanced in the overlaid
%                    image. Default, 1.
%    'startFrameNo': The starting frame number to be merged.
%    'numFrames' : The number of frames to be processed.
%    'edgeDir'   : Path to detected cell edge. If it is specified, cell edge will be plotted. 
%                  Default, ''.
%    'edgeColor' : Color of plotted cell edge. Default, 'r'.
%
% AUTHOR: Lin Ji, Nov. 26, 2005

edgeDir   = '';
edgeColor = 'r';

[redFirstImgFile redInDir filterIndex] = uigetfile('*.*','Pick the first image for red channel');
if isequal(redFirstImgFile,0) || isequal(redInDir,0)
   disp('User pressed cancel. No file is selected.');
   return;
end

[path,redFirstImgFName,no,ext] = getFilenameBody(redFirstImgFile);
redFirstImgFIndex = str2num(no);

%Default output file name.
outFileName = ['rg_' redFirstImgFName];

[redImgFList,index] = getNamedFiles(redInDir,redFirstImgFName);
selIndex = find(index>=redFirstImgFIndex);
redImgFIndex    = index(selIndex);
redImgFList = redImgFList(selIndex);

numRedImgFiles  = length(redImgFList);

[greenFirstImgFile greenInDir filterIndex] = uigetfile('*.*','Pick the first image for green channel');
if isequal(greenFirstImgFile,0) || isequal(greenInDir,0)
   disp('User pressed cancel. No file is selected.');
   return;
end

[path,greenFirstImgFName,no,ext] = getFilenameBody(greenFirstImgFile);
greenFirstImgFIndex = str2num(no);

[greenImgFList,index] = getNamedFiles(greenInDir,greenFirstImgFName);
selIndex = find(index>=greenFirstImgFIndex);
greenImgFIndex    = index(selIndex);
greenImgFList = greenImgFList(selIndex);

numGreenImgFiles  = length(greenImgFList);

if numRedImgFiles ~= numGreenImgFiles
   error('The number of image files in red and green channels do not match.');
end

parentDir = [redInDir filesep '..'];
outDir = uigetdir(parentDir,'Select output directory');

if isequal(outDir,0)
   disp('User pressed cancel. No output directory is selected.');
   return;
end

if ~isdir([outDir filesep 'TIFF'])
   [success,msg,msgID] = mkdir(outDir,'TIFF');
   if ~success
      error('Trouble making directory.');
   end
end
outTifDir = [outDir filesep 'TIFF'];

%Default
redImgIScale     = [0 1];
greenImgIScale   = [0 1];
greenToRedIRatio = 1;
startFrameNo     = 1;
numFrames        = numRedImgFiles;

if nargin > 0
   for kk = 1:2:nargin
      switch varargin{kk}
         case 'outFileName'
            outFileName = varargin{kk+1};
         case 'redImgIScale'
            redImgIScale = varargin{kk+1};
         case 'greenImgIScale'
            greenImgIScale = varargin{kk+1};
         case 'greenToRedIRatio'
            greenToRedIRatio = varargin{kk+1};
         case 'startFrameNo'
            startFrameNo = varargin{kk+1};
         case 'numFrames'
            numFrames = varargin{kk+1};
         case 'edgeDir'
            edgeDir = varargin{kk+1};
         case 'edgeColor'
            edgeColor = varargin{kk+1};
      end
   end
end

%Make sure the requested number of frames does not exceed the total number of images available.
numFrames = min(numFrames,numRedImgFiles-startFrameNo+1);

endFrameNo = startFrameNo+numFrames-1;
rgAVIFile  = [outDir filesep outFileName num2str(startFrameNo) '_' ...
   num2str(endFrameNo) '.avi'];
%To make an avi movie.
aviobj = avifile(rgAVIFile);

pixel_edge = [];
if isdir(edgeDir)
   pixel_edgeFile = [edgeDir filesep 'pixel_edge.mat'];
   if exist(pixel_edgeFile,'file')
      s = load(pixel_edgeFile);
      pixel_edge =s.pixel_edge;
   end
end

figH = figure; hold off;
for kk = startFrameNo:startFrameNo+numFrames-1
   [path,imgFName,no,ext] = getFilenameBody(redImgFList{kk});
   imgIndexStr = no;

   redImg = double(imread([redInDir filesep redImgFList{kk}]));
   redImg = redImg-min(redImg(:));
   rgImg = zeros([size(redImg) 3]);

   maxRedImgI = max(redImg(:));

   redImgF = redImg;
   redImgF(find(redImg>redImgIScale(2)*maxRedImgI)) = redImgIScale(2)*maxRedImgI;
   redImgF(find(redImg<redImgIScale(1)*maxRedImgI)) = redImgIScale(1)*maxRedImgI;

   maxRedImgI = max(redImgF(:));
   minRedImgI = min(redImgF(:));

   redImgF = (redImgF-minRedImgI)/(maxRedImgI-minRedImgI);

   greenImg = double(imread([greenInDir filesep greenImgFList{kk}]));
   greenImg = greenImg-min(greenImg(:));
   rgImg = zeros([size(greenImg) 3]);

   maxGreenImgI = max(greenImg(:));

   greenImgF = greenImg;
   greenImgF(find(greenImg>greenImgIScale(2)*maxGreenImgI)) = greenImgIScale(2)*maxGreenImgI;
   greenImgF(find(greenImg<greenImgIScale(1)*maxGreenImgI)) = greenImgIScale(1)*maxGreenImgI;

   maxGreenImgI = max(greenImgF(:));
   minGreenImgI = min(greenImgF(:));

   greenImgF = (greenImgF-minGreenImgI)/(maxGreenImgI-minGreenImgI);

   if greenToRedIRatio <= 1
      rgImg(:,:,1) = redImgF;
      rgImg(:,:,2) = greenImgF*greenToRedIRatio;
   else
      rgImg(:,:,1) = redImgF/greenToRedIRatio;
      rgImg(:,:,2) = greenImgF;
   end

   figure(figH); hold off;
   imshow(rgImg,[]); hold on;

   if ~isempty(pixel_edge)
      lineObj = plot(pixel_edge{kk}(:,1),pixel_edge{kk}(:,2),edgeColor);
      set(lineObj,'LineWidth',2);
   end
   
   rgImgFile = [outDir filesep outFileName imgIndexStr '.fig'];
   rgImgTifFile = [outTifDir filesep outFileName imgIndexStr '.tif'];
   saveas(figH,rgImgFile,'fig');
   
   %Add frame to avi movie:
   frm = getframe(figH);
   aviobj = addframe(aviobj,frm);
   
   [tifImg,imgCMap] = frame2im(frm);
   imwrite(tifImg,rgImgTifFile,'TIFF','Compression','none');
end

%Close the AVI file.
aviobj = close(aviobj);
