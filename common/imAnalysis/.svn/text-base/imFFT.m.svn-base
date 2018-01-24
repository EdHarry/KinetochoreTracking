function [imgF,varargout] = imFFT(img,varargin)
%imFFT : Time or space Fourier Transform of images.
%
% SYNTAX :
%    imgF = imFFT(img);
%       When 'img' is a cell array of images (>1), do a time domain fourier
%       transform by default. When 'img' is a single image matrix or a cell
%       array of one single image, do a space domain fourier transform by
%       default.
%    imgF = imFFT(img,'FFTAxis','time');
%       Specify a time domain fourier transform. In this case, 'img' has to be
%       a cell array of >1 images.
%    [imgF,f] = imFFT(img,'FFTAxis','time','TimeStep',10);
%       You can also specify a time step (in seconds) and out put the
%       frequency in 'f'.
%
% INPUT :
%    img : Specify images to be analyzed. It can be one single image matrix,
%          a cell array of images or image file names, or a filename of the 
%          first image in a directory. When it is [], a gui window will pop up
%          to select the first image.
% Optional parameters in property/value pairs:
%    Properties :        Values
%    ---------------------------------------------------------------
%    'FFTAxis'  : 'time' or 'space'.
%                 Specify the axis of fourier transform.
%    'TimeStep' : a positive numerical value. Unit: second.
%
% OUTPUT :
%    imgF : Fourier transform of the input image. We only show the positive
%            frequency part.
%    f    : Positive frequencies for time domain.

%Check inputs.
if mod(nargin,2) == 0
   error('The optional parameter/value has to appear in pair.');
end

if nargout > 2
   error('Too many output arguments.');
end

if isempty(img)
   [fName,dirName] = uigetfile('*.tif','imFFT ...');
   img = [dirName,filesep,fName];
end

if ischar(img)
   img = getFileStackNames(img);
end

if iscell(img)
   numImages = length(img);

   %Stack the cell array of images into a 3D matrix with the 3rd dimention 
   % for time.
   if ischar(img{1})
      firstImg = imread(img{1});
      img3D = zeros(size(firstImg,1),size(firstImg,2),numImages);
      for k = 1:numImages
         img3D(:,:,k) = double(imread(img{k}));
      end
   else
      img3D = zeros(size(img{1},1),size(img{1},2),numImages);
      for k = 1:numImages
         img3D(:,:,k) = double(img{k});
      end
   end
else
   numImages = 1;
   img3D(:,:,1) = img;
end

%Scale the image intensity to be from 0 to 1.
maxI = max(img3D(:));
minI = min(img3D(:));

if maxI - minI ~= 0
   img3D = (img3D-minI)/(maxI-minI);
end

if numImages == 1
   FFTAxis   = 'space';
else
   FFTAxis  = 'time';
   TimeStep = 1; %Unit: second.
end

%Check and pass optional parameters.
par = [];
j   = 1;
for k = 1:2:nargin-1
   par{j}   = varargin{k};
   value{j} = varargin{k+1};
   j = j+1;
end

for j = 1:length(par)
   switch par{j}
      case 'FFTAxis'
         if ~ischar(value{j}) | (strcmp(value{j},'time') == 0 & ...
            strcmp(value{j},'space') == 0)
            error('The value for parameter ''interp'' is not valid.');
         end
         FFTAxis = value{j};
      case 'TimeStep'
         if ~isnumeric(value{j}) | value{j} <= 0
            error('The value for parameter ''TimeStep'' is not valid.');
         end
         TimeStep = value{j};
      otherwise
         error(['Parameter' par{j} 'is not recogonized.']);
   end
end

if strcmp(FFTAxis,'time') == 1
   if numImages < 2
      error(['To do time domain fourier transform, there has to be more ' ...
         'than one image.']);
   end

   %Remove DC part.
   imgAvg = sum(img3D,3)/numImages;
   for k = 1:numImages
      img3D(:,:,k) = img3D(:,:,k) - imgAvg;
   end

   %We make the number of sampling points (time steps or number of images in 
   % this case) to be a power of 2 for fast fourier transform. 
   nt2 = 2^ceil(log2(numImages));
   df  = 1/nt2/TimeStep; %Base frequency.
   f   = df*[0:1:nt2/2];

   imgF = (TimeStep/2/pi)*fft(img3D,nt2,3);

   %We only need the positive part of frequencies.
   imgF = imgF(:,:,1:length(f)); 

   if nargout == 2
      varargout{1} = f;
   end
end
