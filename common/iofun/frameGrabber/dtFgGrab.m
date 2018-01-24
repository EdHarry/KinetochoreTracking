function [img, fRoi] = dtFgGrab(roi,rgb)
%DTFGGRAB grabs an image with the DataTranslation framegrabber
%
% SYNOPSIS  [img, fRoi] = dtFgGrab(roi,rgb)
%
% INPUT     roi : (optional) region of interest (ROI) to be grabbed;
%            The ROI will be clipped in case of inconsitency with 
%            the grabbing capabilities. If no ROI is specified, a
%            frame of full size will be grabbed.
%            roi = [upper_left_corner_1, upper_left_corner_2,
%                   width, height]

%           rgb[optional] = 0 : if framegrabber provides rgb frames 
%                               convert to intensity image
%                           1 : if framegrabber provides rgb frames
%                               color information is retained
% OUTPUT    img: either intensity or RGB image
%           fRoi: actual ROI which has been grabbed
%             roi = [upper_left_corner_1, upper_left_corner_2,
%                    width, height]
%
% NOTE      Function refers to the global variable dtFgIsOpen__
%           and runs only if this is set to TRUE

global dtFgIsOpen__;

if dtFgIsOpen__ == 0
   img = [];
   fRoi = [];
   error('DT framegrabber has not been opened => empty return image');
end;	
% call the mexFunction which prints the device info
if(nargin > 0)
   [status, dtFgIsOpen__, result, fRoi] = mexDTFgHandler('grab',roi);
else
   [status, dtFgIsOpen__, result, fRoi] = mexDTFgHandler('grab');
end;

if ndims(result) == 3
   img = permute(result,[3,2,1]);
   [imgHeight,imgWidth,imgDepth] = size(img);
   if imgDepth > 3
      % dealing with an RGB image
		img(:,:,4:imgDepth) = [];
	   if( (nargin < 1) | ~rgb )
   	   auxImg = img;
         img = rgb2gray(auxImg);
      end;
   end;
end;

return;  