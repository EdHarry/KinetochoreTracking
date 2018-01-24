function varargout = imInterp(varargin)
%imInterp Interpolation of image intensity to subpixel points and return the 
%         intensity as double values.
%
% SYNOPSIS :
%    outI = imInterp(img,YX)
%    outI = imInterp(img,{y x})
%    The function interpolates the intensity of an image to points specified
%    in 'YX' (or 'y' and 'x'). The output of the interpolated values are of 
%    class double and are linearly scaled to values from 0 to 1. 
%
%    outI = imInterp(img)
%    If the interpolation points are not specified the interpolation is 
%    carried on all the pixels in 'img'.
%
%    outI = imInterp(img,[],'Gaussian',...)
%    outI = imInterp(img,YX,method,par1,...,parm)
%    outI = imInterp(img,{y x},method,par1,...,parm)
%    An interpolation method can also be used. See INPUT below for available
%    interpolation methods. By default, the 'Gaussian' convolution method is 
%    used.
%
%    [sp,outI] = imInterp(img,YX,method,par1,...,parm)
%    When 'method' is the spline interpolation is used, the pp- or B-form of 
%    the interpolation, can also be output in the variable, 'sp'. 
%
%    outI = imInterp(img,YX,'spap2',...)
%    sp   = imInterp(img,[],'spap2',...)
%    When no interpolation point is specified and 'method' is one of the 
%    spline interpolation methods, the output is the B-form or pp-form of
%    the spline interpolation. When the interpolation points are specified,
%    the output is the value itself.
%
% INPUT :
%    img : The image intensity matrix. It can be of class uint8, uint16, or
%       double. It can also be a string that specifies the name of the image
%       file.
%    YX  : An m-by-2 matrix that specifies the coordinates of m points where
%       we want to find the values by interpolation. It is given in the form
%       [y x] where 'x' is the x coordinates and 'y' is the y coordinates. The
%       unit of the coordinates are in pixels. For points that are outside of
%       the image, the output is NaN.
%    {y x} : The interpolated points in this case are the grid points generated
%       by '[Y,X] = ndgrid(y,x)'.
%    method : A string that specifies the interpolation method to use.
%       'spap2'   : Least square B-Spline.
%       'spaps'   : Smoothing spline.
%       'Gausian' : Convolution with Gausian kernel.
%    par1,...,parm : Parameters that are needed for the various interpolation
%       methods.
%        method        (par1,...,parm)
%       ------------------------------------
%       'spap2'   : ({knotsY,knotsX}, [kY kX])
%                   ([dY dX], [kY kX])
%                   (dd, k) %In this case, dY=dX=dd, kY=kX=k.
%       'spaps;   : (tol)
%       'Gaussian' : (corLen) %In this case, corLenY=corLenX=corLen.
%                   ([corLenY corLenX])
%
% OUTPUT :
%    outI : The scaled intensities at the points given in 'YX' (or 'y' and
%       'x'). It is of size m-by-1 in the case of 'YX' (m-by-2 matrix) and of
%       size length(y)-by-length(x) in the case of 'y' and 'x'. The values are
%       of class double even if the input image is of class uint8, or uint16.
%    sp : The B-form or pp-form of the spline interpolation when 'method' is 
%       one of the spline interpolation methods.
%
% AUTHOR : Lin Ji 
% DATE   : Feb. 12, 2004.
% Changed by Kathryn to convert pixel xy-coordinates correctly
% outI = outI((floor(YX(:,2)+0.5)-1)*numPixelsY+floor(YX(:,1)+0.5)); 

if nargin < 1
   error('There has to be at least one input argument.');
end

img = varargin{1};

%Check the legitimacy of the input, 'img'.
if ischar(img)
   imI = imread(img);
else
   imI = img;
end

if size(imI,3) == 3
   error('RGB images are not supported. Call RGB2GRAY first.');
%elseif isind(imI)
%   error('Indexed images are not supported. Call IND2GRAY first.');
%elseif ~isgray(imI) | ~isind(imI)
%   error('Not an image intensity matrix.');
end

numPixelsX = size(imI,2);
numPixelsY = size(imI,1);

if nargin >= 3 & ~isempty(varargin{3})
   method = varargin{3};
else
   %Default interpolation method.
   method = 'Gaussian';
end

YX = []; y = []; x = [];
if nargin >= 2 
   if ~isempty(varargin{2})
      if iscell(varargin{2})
         y = varargin{2}{1};
         x = varargin{2}{2};
         %Check if the interpolation points are in the ranage of the image.
         if min(y) < 1 | max(y) > numPixelsY | min(x) < 1 | max(x) > numPixelsX
            error('The interpolation points are out of the image range.');
         end
      else
         YX = varargin{2}; 
      end
   end
end

%Check if 'YX' or 'y' and 'x' are correctly defined.
if ~isempty(YX) & (~isnumeric(YX) | ndims(YX) > 2 | size(YX,2) ~= 2)
   error(['The interpolation points is not correctly defined. ' ...
      'See help imInterp.']);
elseif (~isempty(x) & isempty(y)) | (~isempty(y) & isempty(x)) 
   error(['The interpolation points is not correctly defined. ' ...
      'See help imInterp.']);
elseif ~isnumeric(y) | ~isnumeric(x) | ...
   ndims(y) > 2 | ndims(x) > 2 | min(size(y)) > 1 | min(size(x)) >1
   error(['The interpolation points is not correctly defined. ' ...
      'See help imInterp.']);
end

%Parse the rest of the arguments according to 'method'.
if strcmp(method,'Gaussian') == 1
   if nargout > 1
      error('Too many output arguments.');
   end

   %Default correlation length.
   corLen = 0.1; 

   if nargin > 4
      error(['Too many input arguments for the specified ' ...
         'interpolation method.']);
   end

   if nargin == 4 & ~isempty(varargin{4})
      corLen = varargin{4};
   end

   %Check if 'corLen' is set correctly.
   if ~isnumeric(corLen) | (isnumeric(corLen) & length(corLen) > 2)
      error(['The correlation length should be provided as a ' ...
         'numerical array of 1 or 2 numbers.']);
   elseif min(corLen) < 0
      error('The correlation length has to be positive.');
   end

   if length(corLen) == 1
      corLenY = corLen;
      corLenX = corLen;
   else
      corLenY = corLen(1);
      corLenX = corLen(2);
   end
elseif strcmp(method,'spap2') == 1
   if nargout > 2
      error('Too many output arguments.');
   end

   order = 4; %Default order of spline interpolation.

   if nargin > 5 
      error(['Too many input arguments for the specified ' ...
         'interpolation method.']);
   end

   if nargin == 5 & ~isempty(varargin{5})
      order = varargin{5};
   end

   %Check if 'order' is correctly defined.
   if ~isnumeric(order) | (isnumeric(order) & length(order) > 2)
      error(['The order of the spline interpolation should be ' ...
         'provided as a numerical array of 1 or 2 positive integers.']);
   elseif min(order) < 0
      error('The order of the spline interpolation has to be positive.');
   end

   if length(order) == 1
      orderY = order;
      orderX = order;
   else
      orderY = order(1);
      orderX = order(2);
   end

   if floor(orderX) ~= orderX | floor(orderY) ~= orderY
      error('The order of the spline interpolation has to be integer.');
   end

   if nargin >= 4 & ~isempty(varargin{4})
      if iscell(varargin{4}) & length(varargin{4}) == 2
         knotsY = varargin{4}{1};
         knotsX = varargin{4}{2};
      elseif isnumeric(varargin{4}) 
         %Get distance between breaks for constructing knot sequence.
         if length(varargin{4}) == 2
            dY = varargin{4}(1);
            dX = varargin{4}(2);
         elseif length(varargin{4}) == 1
            dY = varargin{4};
            dX = dY;
         else
            error(['The arguments for the method, ''spap2'' ' ...
               'are not correctly given.']);
         end

         numBrksY = max(3,ceil(numPixelsY/dY));
         numBrksX = max(3,ceil(numPixelsX/dX));
         knotsY   = augknt(linspace(1,numPixelsY,numBrksY),orderY);
         knotsX   = augknt(linspace(1,numPixelsX,numBrksX),orderX);
      else
         error(['The arguments for the method, ''spap2'' ' ...
            'are not correctly given.']);
      end
   else
      %If no knot sequence is provided, the default knot sequence with
      % intervals of appoximately every other pixel is created and the default
      % order is 4.
      numBrksY = max(3,ceil(numPixelsY/2));
      numBrksX = max(3,ceil(numPixelsX/2));
      knotsY   = augknt(linspace(1,numPixelsY,numBrksY),orderY);
      knotsX   = augknt(linspace(1,numPixelsX,numBrksX),orderX);
   end
elseif strcmp(method,'spaps') == 1
   if nargout > 2
      error('Too many output arguments.');
   end

   %Default correlation length.
   tol = 0.01; 

   if nargin > 4
      error(['Too many input arguments for the specified ' ...
         'interpolation method.']);
   end

   if nargin == 4 & ~isempty(varargin{4})
      tol = varargin{4};
   end

   %Check if 'corLen' is set correctly.
   if ~isnumeric(tol) | length(tol) > 1
      error(['The tolerance for smoothing spline should be provided as a ' ...
         'single numerical value.']);
   elseif tol < 0
      error('The tolerance has to be positive.');
   end
else
   error('The specified interpolation method is not recogonized.');
end

%Rescale the image intensity to the range [0 1]. First, convert 
% 'imI' to double if it is of class 'uint8' or 'uint16'.
imI  = double(imI);
minI = min(imI(:));
maxI = max(imI(:));
imI  = (imI-minI)/(maxI-minI);

%Start the interpolation.
if strcmp(method,'Gaussian') == 1
   %Compute the matrix of weights for intensities in a correlation square
   % whose center is moved around to each interpolation point. The size of 
   % the sqare is determined by 'corLenY' and 'corLenX'.
   corY = [-floor(2*corLenY):floor(2*corLenY)].';
   corX = [-floor(2*corLenX):floor(2*corLenX)];
   if length(corY) == 1 & length(corX) == 1
      %Exact interpolation in this case.
      outI = imI;
   else
      W    = exp(-corY.^2/2/corLenY^2)*exp(-corX.^2/2/corLenX^2);
      W    = W/sum(sum(W));

      outI = conv2(imI,W,'same');
   end

   if ~isempty(YX)
      outI = outI((floor(YX(:,2)+0.5)-1)*numPixelsY+floor(YX(:,1)+0.5));

      %The following for loop should be implemented as mex function. The
      % initial idea is that when we only ask for the intensity values at a few
      % points, we don't need to convolve the whole image. Anyone has an idea 
      % to avoid the for loop here or is interested in implementing the mex
      % function is highly respected.

      %outI = zeros(size(YX));
      %for k = 1:size(YX,1)
      %   %Move the center of the correlation square defined by 'corY' and 
      %   % 'corX' to each interpolation point.
      %   corYM = corY+floor(YX(k,1)+0.5);
      %   corXM = corX+floor(YX(k,2)+0.5);

      %   %The square 'corYM x corXM' has to be inside the image.
      %   indY = find(corYM>0&corYM<=numPixelsY);
      %   indX = find(corXM>0&corXM<=numPixelsX);
      %   outI(k) = sum(sum(imI(corYM(indY),corXM(indX)).*W(indY,indX)))/ ...
      %      sum(sum(W(indY,indX)));
      %end
   elseif ~isempty(x)
      outI = outI(floor(y+0.5),floor(x+0.5));
      %outI = zeros(length(y),length(x));
      %for jy = 1:0 %length(y)
      %   for jx = 1:length(x)
      %      %Move the center of the correlation square 'corY x corX'
      %      % to each interpolation point.
      %      corYM = corY+floor(y(jy)+0.5);
      %      corXM = corX+floor(x(jx)+0.5);

      %      %The square 'corYM x corXM' has to be inside the image.
      %      indY = find(corYM>0&corYM<=numPixelsY);
      %      indX = find(corXM>0&corXM<=numPixelsX);
      %      outI(jy,jx) = sum(sum(imI(corYM(indY),corXM(indX)).* ...
      %         W(indY,indX)))/sum(sum(W(indY,indX)));
      %   end
      %end
   end
   varargout{1} = outI;
   return;
end

if strcmp(method,'spap2') == 1
   sp = spap2({knotsY,knotsX},[orderY orderX], ...
      {[1:numPixelsY],[1:numPixelsX]},imI);
elseif strcmp(method,'spaps') == 1
   sp = spaps({[1:numPixelsY],[1:numPixelsX]},imI,tol);
end

if nargout == 2
   varargout{1} = sp;
   if ~isempty(YX)
      varargout{2} = fnval(sp,YX.');
   elseif ~isempty(x)
      varargout{2} = fnval(sp,{y,x});
   else
      varargout{2} = fnval(sp,{[1:numPixelsY],[1:numPixelsX]});
   end
else
   if isempty(YX) & isempty(x)
      varargout{1} = sp;
   elseif ~isempty(YX)
      varargout{1} = fnval(sp,YX.');
   elseif ~isempty(x)
      varargout{1} = fnval(sp,{y,x});
   end

end

