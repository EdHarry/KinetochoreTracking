function VPanelResize(figHandle,offset)
%VPANELRESIZE sets axes according to window size
%
% SYNOPSIS VPanelResize(figHandle, offset)
%
% INPUT figHandle: a figure handle
%		  offset   : (optional) a rect which gives the offset from [top bottom left right]
% c: 14/6/99	dT

if(nargin==0)
   figH=gcbf;
   offset=[0 0 0 0];  
elseif (nargin==1)
   figH = figHandle;
   offset=[0 0 0 0];
else
   figH = figHandle;
end;

axesH=findobj(figH,'type','axes');
%no resize if multiples axes
if length(axesH)>1
    %axis(axesH,'image');
    return;
end;

axesH = get(figH,'CurrentAxes');
%set to units to pixels!
set(axesH,'units','pixels');
panSize = get(figH,'Position');
if (((panSize(4)-(offset(1)+offset(2)))>0) & ((panSize(3)-(offset(3)+offset(4)))>0))
   set(axesH,'Position',[offset(3) offset(2) (panSize(3)-(offset(3)+offset(4)))...
         (panSize(4)-(offset(1)+offset(2)))]);
end;
