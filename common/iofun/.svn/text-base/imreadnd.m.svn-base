function [I,adepth,msg,res]=imreadnd(filename,adepth,algorithm)
%IMREADNB reads an image file and normalizes all pixel intensities 
%   to the range 0..1 (the information is stored in matlab 'double' 
%   precision type).
%
%   SYNOPSIS [I,adepth,msg]=IMREADND(filename,adepth,algorithm)
%
%   INPUT   filename:   name of the file with entire path
%           adepth:     number of bits actually used to represent 
%                       image information (enter '-1' for automatic 
%                       search)
%
%           algorithm:  0 = 'default' 
%                       1 = 'smart' 
%
%                       This parameter will be considered only if adepth 
%                       has been set to '-1'
%
%   OUTPUT  I:          double image matrix normalized to 0..1
%           adepth:     the actual depth given by the user or calculated
%                       is returned
%           msg:        info/error string returned by imreadnd 
%           res:        0: success 
%                       -1: error
%
%   NOTE: if only the parameter 'filename' is given, the function will
%         calculate adepth using the 'default' algorithm
%
%
%
%   IMREADNB relies on matlab function IMREAD, which can load
%   following data types:
%
%      Extension       File type
%   --------------------------------------------------------------
%      'jpg' or 'jpeg' Joint Photographic Experts Group (JPEG)
%      'tif' or 'tiff' Tagged Image File Format (TIFF)
%      'bmp'           Windows Bitmap (BMP)
%      'png'           Portable Network Graphics
%      'hdf'           Hierarchical Data Format (HDF)
%      'pcx'           Windows Paintbrush (PCX)
%      'xwd'           X Window Dump (XWD)
%
%   Type HELP IMREAD to get more accurate information about file
%   formats that can be read.
%
%   Use IMREAD instead if you only want to load an image into 
%   memory without modifying its content. 
%   And remember that IMREAD does not take the variable 'adepth'
%   as a parameter.
%

if nargin==0
   error('This functions requires at least one parameter');   
end
if nargin==1
   alg=[1 0];
end
if nargin==2
   if adepth==-1
      alg=[1 0];
   else
      alg=[0 0];
   end
end
if nargin==3
   if adepth==-1
      alg=[1 algorithm];
   else
      alg=[0 0];
   end
end

% 'filename' is passed 'as is' to imread, which will check it
[X] = imread(filename);

% Conversion of I to type double (necessary for normalization)
X=double(X);

% *************************
%
% EVALUATION
%
% *************************

% If adepth has been defined, normalization is performed to the given depth
if alg==[0 0]
   % Setting xmin and xmax
   xmin=0;
   xmax=2^adepth-1;
end
if alg==[1 0]
   % 'Default' algorithm: looking for maximun intensity
   tmax=max(max(X));
   % Setting xmin
   xmin=0;
   % Calculating xmax for normalization
   if tmax~=0
      if (log2(tmax)-fix(log2(tmax)))>0
         adepth=fix(log2(tmax))+1;
      else
         adepth=log2(tmax);   
      end
      xmax=2^adepth-1;
   else
      xmax=0;
   end
end
if alg==[1 1]
   % 'Smart' algorithm: looking for minimum and maximun intensity
   tmin=min(min(X));
   tmax=max(max(X));
   % Calculating xmin for normalization
   if tmin~=0
      xmin=2^(fix(log2(tmin)));
   else
      xmin=0;
   end
   % Calculating xmax for normalization   
   if tmax~=0
      if (log2(tmax)-fix(log2(tmax)))>0
         adepth=fix(log2(tmax))+1;
      else
         adepth=log2(tmax);
      end
      xmax=2^adepth-1;
   else
      xmax=0;
   end
end

% Normalization
I=(X-xmin)/(xmax-xmin);

% Check for normalization
ch=max(max(I));
if ch>1
   msg='ERROR: Normalization failed! Select a greater bit depth.';
   res=-1;
else
   msg=strcat('The image has been normalized to:',num2str(adepth),' bit.');
   res=0;
   if alg(2)~=-1
   end
end

