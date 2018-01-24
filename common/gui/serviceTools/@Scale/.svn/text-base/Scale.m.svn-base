function sc=Scale(handle)
%SCALE create a scale in a image
%
% SYNOPSIS sc=Scale(handle)
%
% INPUT  handle: an image handle OR a figure handle
%
% constructor:
%                   Scale(handle)
%
% public methods:
%                   SetScale(scale)
%
% private methods:
%                   OptionWindow(scale)	
%
% c: 12/8/99	dT

% Check whether we got an image or a figure and initialize scale
if (strcmp(get(handle,'Type'),'image'))
   imAxesH = get(handle,'Parent');
   sc.figure = get(imAxesH,'Parent');
else
   sc.figure = handle;
end;
%scale length in mu
sc.length = 100;
% origin of the scale in the figure
sc.origin=[];
% one pixel corresponds to ... mu
sc.pixelSize =2;
% line Handle
sc.lineH = [];
% text Handle
sc.textH = [];
% additional text
sc.text=[num2str(sc.length) ' \mum'];
% color (1 = white, 2 = black)
sc.color = 1;
sc = class(sc,'Scale');
