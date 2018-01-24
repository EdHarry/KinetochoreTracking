function imT=ImAdjustTool(handle)
%IMADJUSTTOOL gamma correction tool for images
%
% SYNOPSIS imT=ImAdjustTool(handle)
%
% INPUT  handle: an image handle OR a figure handle
%
% constructor:
%                   ImAdjustTool(handle)
%
% public methods:
%                   DoImageAdjust
%
% private methods:
%                   ImageAdjustToolWindow
%
% c: 17/8/99	dT

% Check whether we got an image or a figure and initialize gammaobj
if (strcmp(get(handle,'Type'),'image'))
   imAxesH = get(handle,'Parent');
   imT.figure = get(imAxesH,'Parent');
else
   imT.figure = handle;
end;

% Check if there's an imAdjustTool window open which is pointing to this figure
auxH=findobj(0,'Tag','IMAGE_ADJUST_TOOL');
if (~isempty(auxH))
   for l = 1:length(auxH)
   	auxImT=get(auxH(l),'UserData');
      if(imT.figure==auxImT.figure)
         %we found a valid imAdjustTool
         imT=auxImT;
         return;
      end;
   end;
end;

% initialize ImTool
figure(imT.figure);
imT.oldCMap = colormap;
imT.newCMap=[];
imT.oldInt=[0 1];
imT.newInt=[0 1];
imT.gamma=1;
% handle to the tool panel
imT.panelH=[];
   
imT = class(imT,'ImAdjustTool');
   


