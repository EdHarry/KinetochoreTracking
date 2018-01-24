function DoImageAdjust(imT)
%IMADJUSTTOOL/DOIMAGEADJUST do adjust for object 'imT'
%
% SYNOPSIS DoImAdjust(imT)
%
% INPUT  imT: an imageToolAdjust object
%
% c: 18/8/99	dT

if (isempty(imT.panelH))
   % Open the imAdjustTool window
   imT.panelH=ImageAdjustToolWindow(imT);   
   % Connect the Scale to the options window   
   set(imT.panelH,'UserData',imT);
else
   figure(imT.panelH);
end;
