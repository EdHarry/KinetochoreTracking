function [img,fRoi] = grab(roi)
%GRAB grabs a frame with the globally specified framegrabber board
%
% The type of framegrabber is read from the global variable fgType__
%
% SYNOPSIS [img, fRoi] = grab(roi)
%
% INPUT roi: (optional) region of interest (ROI) to be grabbed;
%            The ROI will be clipped in case of inconsitency with 
%            the grabbing capabilities. If no ROI is specified, a
%            frame of full size will be grabbed.
%            roi = [upper_left_corner_1, upper_left_corner_2,
%                   width, height]
%
% OUTPUT img : acquired frame
%        fRoi: actual ROI which has been grabbed
%              roi = [upper_left_corner_1, upper_left_corner_2,
%                     width, height]
%

global fgType__;

% set the grab command
if(strcmp(fgType__,'DT'))
   if(nargin == 1)	
      grabCmd = 'dtFgGrab(roi)';
   else
      grabCmd = 'dtFgGrab';
   end;
end;

% grab the image
[img,fRoi] = eval(grabCmd);
