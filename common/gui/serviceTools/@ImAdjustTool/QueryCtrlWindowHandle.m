function handle = QueryCtrlWindowHandle(imT)
%IMADJUSTTOOL/QUERYCTRLWINDOWHANDLE gets handle to ImAdjustTool control window
%
% SYNOPSIS QueryCtrlWindowHandle(imT)
%
% INPUT  imT: an ImAdjustTool object
% 
% OUTPUT handle : figure handle to the window;
%                 [] if there is none
%
% c: 29/8/00	 gD

handle = imT.panelH;
