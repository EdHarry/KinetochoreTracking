function varargout = imStackTDFilter(stkIn,freqThreshold,varargin)
%imStackTDFilter : Filter a stack of images (stacked as a time series ) in the 
%                  time domain.
%
% SYNTAX :
%    stkOut = imStackTDFilter(stkIn,freqThreshold);
%    stkOut = imStackTDFilter(stkIn,freqThreshold,'outFile',fileName);
%    imStackTDFilter(stkIn,freqThreshold,'outFile',fileName);
%
% INPUT :
%    stkIn : The stack of images to be filtered. It can be a cell array of
%            images, or image file names, or the file name of the first image
%            in a directory. When it is [], a window will pop up to select the
%            first image.
%    freqThreshold : Specify the cut-off frequency in term of the ratio of
%            frequencies (>0 & <1) to be kept in the whole range.
%
% Optional parameters in property/value pairs:
%    Properties :        Values
%    ---------------------------------------------------------------
%    'outFile' : If a common name of the output files is not given, a window
%                will pop up for the output directory.
%
% OUTPUT :
%    stkOut : A cell array of filtered images. When it is not given, a
%             directory has to be specified to save the filtered images.

if nargout > 1
   error('Too many output arguments. See help imStackTDFilter.');
end

%Check inputs.
if mod(nargin,2) ~= 0
   error('The optional parameter/value has to appear in pair.');
end

if freqThreshold <= 0 | freqThreshold >= 1
   error(['The value of ''freqThreshold'' that specifies the proportion ' ...
      '''of frequencies to be kept is between 0 and 1.']);
end
if isempty(stkIn)
   [fName,dirName] = uigetfile('*.tif','Select the first image file ...');
   stkIn = [dirName,filesep,fName];
end

if ischar(stkIn)
   stkIn = getFileStackNames(stkIn);
end

numImages = length(stkIn);

%Stack the cell array of images into a 3D matrix with the 3rd dimention 
% for time.
if ischar(stkIn{1})
   firstImg = imread(stkIn{1});
   img3D = zeros(size(firstImg,1),size(firstImg,2),numImages);
   for k = 1:numImages
      img3D(:,:,k) = double(imread(stkIn{k}));
   end
else
   img3D = zeros(size(stkIn{1},1),size(stkIn{1},2),numImages);
   for k = 1:numImages
      img3D(:,:,k) = double(stkIn{k});
   end
end

%Default optional parameters.
outFile = [];

%Check and pass optional parameters.
par = [];
j   = 1;
for k = 1:2:nargin-2
   par{j}   = varargin{k};
   value{j} = varargin{k+1};
   j = j+1;
end

for j = 1:length(par)
   switch par{j}
      case 'outFile'
         if ~ischar(value{j})
            error('The value for parameter ''outFile'' is not valid.');
         end
         outFile = value{j};
      otherwise
         error(['Parameter' par{j} 'is not recogonized.']);
   end
end

%Fourier transform the data images in time.
if numImages < 2
   error(['There have to be at least two images in the stack.']);
end

%Scale the image value range to [0 1];
minI  = min(img3D(:));
maxI  = max(img3D(:));
img3D = (img3D-minI)/(maxI-minI);

nt2 = 2^ceil(log2(numImages));
imgF = fft(img3D,nt2,3);

%Filtering by cutting off high frequencies.
cutOffIndex = ceil(nt2/2*freqThreshold);
imgF(:,:,cutOffIndex+1:nt2-cutOffIndex+1) = 0; 

img3D = real(ifft(imgF,nt2,3));

%Scale the image value range to [0 1];
minI  = min(img3D(:));
maxI  = max(img3D(:));
img3D = (img3D-minI)/(maxI-minI);

%Convert it to uint16 image.
img3D = uint16(round(img3D*65535));

if nargout == 1
   for k = 1:numImages
      varargout{1}{k} = img3D(:,:,k);
   end
end

if isempty(outFile)
   return;
end

% Select a directory to save filtered images.
[file,outDir]=uiputfile('Click on SAVE','Select output directory');
if outDir==0
    disp('No image files are saved.');
    return
end    

fprintf(1,'Saving ...');
indexFormat = sprintf('%%.%dd',length(num2str(numImages)));
for k = 1:numImages
   indexStr = sprintf(indexFormat,k);
   imwrite(img3D(:,:,k),[outDir,filesep,outFile,indexStr,'.tif']);
end

