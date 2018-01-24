function [stack,fRoi] = grabSeq(nImg,freq,roi)
%GRABSEQ grabs a sequence
%
% This function stores the stack in RAM. Thus, only a limited
% number of images should be grabbed. The type of framegrabber 
% is read from the global variable fgType__
%
% SYNOPSIS [stack,fRoi] = grabSeq(nImg,freq,roi)
%
% INPUT nImg : number of images to be grabbed
%       freq : frequency in [Hz]
%       roi  : (optional) region of interest (ROI) to be grabbed;
%              The ROI will be clipped in case of inconsitency with 
%              the grabbing capabilities. If no ROI is specified, a
%              frame of full size will be grabbed.
%              roi = [upper_left_corner_1, upper_left_corner_2,
%                     width, height]
%
% OUTPUT stack : acquired sequence
%        fRoi  : actual ROI which has been grabbed
%                roi = [upper_left_corner_1, upper_left_corner_2,
%                       width, height]


global fgType__;

% set the grab command
if(strcmp(fgType__,'DT'))
   if(nargin == 3)
      grabCmd = 'dtFgGrab(roi)';
   else	
      grabCmd = 'dtFgGrab';
   end;
else
   %use left camera only
   grabCmd = 'ppMonoFgGrab(0)';
end;

% get time between frames
dTime = 1/abs(freq);

% make sure that no other waitbar is active in the same root,
% this might cause interferences with strange errors
% GD-6-24-98
hWaitbar = findobj(0,'Type','figure','Tag','TMWWaitbar');
if(~isempty(hWaitbar))
   close(hWaitbar);
end;

% start the waitbar
hWaitBar = waitbar(0,'Acquiring ...');
t0 = clock;

% grab the first image
[stack, fRoi] = eval(grabCmd);
waitbar(1/nImg);

for iImg = 2:nImg
   nEtime = (iImg - 1)*dTime;
   % waiting loop
   while(etime(clock,t0) < nEtime)
   end;
   %%%% DEBUG command: etime(clock,t0)
   stack = cat(3,stack,eval(grabCmd));
   % display image if required
   waitbar(iImg/nImg);
end;
close(hWaitBar);
