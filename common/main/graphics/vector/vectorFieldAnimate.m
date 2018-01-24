function M = vectorFieldAnimate(varargin)
%vectorFieldAnimate  Show the animation of the given vector field over time.
%
% SYNOPSIS :
%    M = vectorFieldAnimate(XY,V,scale,'name1',value1,...)
%    M = vectorFieldAnimate(XY,V,...)
%
% INPUT :
%    XY : Points where the vector field is shown. It is an m-by-2 matrix in
%       the form [x y] where 'm' is the number of points.
%    V : The stack of vector fields over time. It is an m-by-2-by-n
%       multidimentional array where 'm' equals the number of points stored in
%       'XY' and 'n' is the number of time steps (or frames). The two colums
%       of the 2nd dimention has the form [vx vy] where 'vx' and 'vy' are the
%       x and y components of the vector respectively.
%    scale : A non-negative numerical value. When it is zero, the vectors
%       will be automatically scaled. The scaling is the same for all the
%       frames in order to display the dynamical change of the vector length.
%       The default value is zero (automatical scaling).
%
%    The following properties can also be specified as optional inputs.
%    'bgImg' : It can be either one image or a cell array of
%       images whose length equals the number of time steps.
%    'aviMovie' : A string or a cell array with 2 elements. In the case of a
%       cell array, the first element is a string that specifies the name of 
%       an AVI movie for output and the second element is itself a cell array 
%       of parameters to be passed to MATLAB function 'movie2avi'. If there is
%       no parameter, a string can be passed for the name of the AVI movie.
%       See help movie2avi.
%    'imgOutDir' : A string that pecifies a directory where the image of each 
%       frame can be saved so that a movie (e.g. QuickTime movie) can be made 
%       afterwards. (Not yet available)
%    'vColor' : The color of the vector. It can be one of the following 
%       values:
%          Color names : such as 'r', 'g' or 'cyan' etc.
%          RGB-triple  : [0.2 0.3 0].
%          Cell array  : of colors in one of the above forms.
%       If only one color is specified, all the fields are drawn in the same
%       color. If a cell array is given, each color in the array is used for 
%       vector field and the colors are cycled when there are less colors
%       than the number of fields.
%    'colorMap' : Either an m-by-3 matrix whose values are between 0 and 1 or
%       one of the matlab accepted color map names such as 'default', 'jet' or
%       'cool' etc. When a color map is specified, each field will be assigned
%       a color that is determined by an even distribution of the vector field
%       indices over the range of the color map. See help colormap.
%    If neither 'vColor' nor 'colorMap' is specified, the default color map
%    'colormap('default')' will be used. If both are specified, 'colorMap' has
%    priority.
%
% OUTPUT :
%    M : A matlab movie is returned.

%Parse the inputs.
if nargin < 2
   error('Not enought input arguments.');
end

XY = varargin{1};
V  = varargin{2};

%Check the legitimacy of 'XY' and 'V'
if ~isnumeric(XY) | ndims(XY) > 2 | size(XY,2) ~= 2
   error(['The first argument should be an m-by-2 numerical matrix that ' ...
      'specifies the coordinates of the points where the vector field ' ...
      'is shown.']);
end

if ~isnumeric(V) | ndims(V) > 3 | size(V,2) ~= 2 | ...
   size(V,1) ~= size(XY,1)
   error(['The second argument that specifies the vector field is not ' ...
      'defined correctly.']);
end

numFrames = size(V,3);

scale  = 0; %Default scale.
pStart = 3; %Where the properties start.
if nargin >= 3 & isnumeric(varargin{3})
   scale = varargin{3};

   if length(scale) ~= 1
      error(['The scale specified in the third argument is ' ...
         'a single numerical value.']);
   end
   if scale < 0
      error(['The scale specified in the third argument can not ' ...
         'be negative.']);
   end
   pStart = 4;
end

%Default property values.
bgImg        = [];
aviMovieName = [];
movie2aviPar = [];
imgOutDir    = [];
vColor       = [];
colorMap     = [];

if nargin >= pStart
   if rem(nargin-pStart+1,2) ~= 0
      error('The properties and values specified do not match as pairs.');
   end

   numProperties = (nargin-pStart+1)/2;
   inProperty    = cell(numProperties,1);
   inValues      = cell(numProperties,1);
   for k = 1:numProperties
      inProperty{k} = varargin{pStart+2*k-2};
      inValue{k}    = varargin{pStart+2*k-1};

      switch inProperty{k}
      case 'bgImg'
         if iscell(inValue{k})
            if length(inValue{k}) ~= 1 | ...
               length(inValue{k}) ~= numFrames
               error(['The number of background images can either be ' ...
                  'one or equal the number of frames.']);
            end

            if ischar(inValue{k}{1})
               for jj = 1:length(inValue{k})
                  bgImg{jj} = imread(inValue{k}{jj});
               end
            else
               bgImg = inValue{k};
            end
         else
            if ischar(inValue{k})
               bgImg{1} = imread(inValue{k});
            else
               bgImg{1} = inValue{k};
            end
         end
      case 'aviMovie'
         if ischar(inValue{k})
            aviMovieName = inValue{k};
         elseif ~iscell(inValue{k}) | length(inValue{k}) ~= 2 | ...
            ~ischar(inValue{k}{1}) | ~iscell(inValue{k}{2})
            error(['The value for Property ''aviMovie'' is not ' ...
               'correctly defined.']);
         else
            aviMovieName = inValue{k}{1};
            movie2aviPar = inValue{k}{2};
         end
      case 'imgOutDir'
         imgOutDir = inValue{k};
      case 'vColor'
         vColor = inValue{k};
      case 'colorMap'
         if ischar(inValue{k})
            %We do it in this complicated way becaus of the way 'colormap'
            % works in matlab.
            %Get the current map.
            curMap   = colormap;
            %Switch to the specified map and assign it to 'colorMap'.
            colormap(inValue{k});
            colorMap = colormap;
            %Restore the original color map.
            colormap(curMap);
         else
            colorMap = inValue{k};
         end
      otherwise
         error([inProperty{k} ' is not an acceptable property.']);
      end
   end
end

if isempty(vColor) & isempty(colorMap)
   curMap = colormap;
   colormap('default');
   colorMap = colormap;
   colormap(curMap);
end

if ~isempty(colorMap)
   %'colorMap' has higher priority.
   %Evenly distributing the indices of the vector fields over the range of the
   % color map. 'cI' is the index vector into the color map.
   cI = floor(linspace(1,size(colorMap,1),numFrames)+0.5);
   %'vC' : the list of colors used for each frame.
   vC = colorMap(cI,:);
else
   if iscell(vColor)
      cLen = length(vColor);

      %When the number of colors, 'cLen' is less than 'numFrames', the
      % colors are reused periodically.
      vC(1:numFrames) = vColor(rem([0:numFrames-1],cLen)+1);
   else
      vC(1:numFrames) = {vColor};
   end
end

%Display the vector field and make the movie.
figure(gcf); hold off;
for k = 1:numFrames
   if ~isempty(bgImg)
      if length(bgImg) == 1
         imshow(bgImg{1},[]); hold on;
      else
         imshow(bgImg{k},[]); hold on;
      end
   end

   h = quiver(XY(:,1),XY(:,2),V(:,1,k)*scale,V(:,2,k)*scale,0);
   if iscell(vC)
      set(h(1),'Color',vC{k});
      set(h(2),'Color',vC{k});
   else
      set(h(1),'Color',vC(k,:));
      set(h(2),'Color',vC(k,:));
   end
   title(sprintf('Scale : %5.2f, Number of Frames : %d',scale,numFrames));
   hold off;
   M(k) = getframe;
end

if ~isempty(aviMovieName)
   if isempty(movie2aviPar)
      movie2avi(M,aviMovieName);
   else
      movie2avi(M,aviMovieName,movie2aviPar{:});
   end
end
