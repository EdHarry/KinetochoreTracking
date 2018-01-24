function ans = uiSSDTrackPanelGetDflt(request)
% service function to get some default values for 
% the SSD track panel

% possible requests : 'maxImgInStack' -> max number of images in stack
%                     'loadNImgToStack' -> number of images load to the stack at one time
%                     'logfname' -> name of the log file
%                     'bdFactor' -> additional border for template and 
%                                   search image: bdFactor*sigma

maxImgInStack = 10;
%loadNImgToStack = 10;
logfname = 'ssdTrack.log';
bdFactor = 5;

switch(request)
case 'maxImgInStack',
   ans = maxImgInStack;
case 'loadNImgToStack',
   ans = loadNImgToStack;
case 'logfname',
   ans = logfname;
case 'bdFactor';
   ans = bdFactor;
end;

