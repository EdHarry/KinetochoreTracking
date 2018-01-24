function ImageProfile(imgH)
%IMAGEPROFILE creates the color profile along user selected sections
% and shows it in a seperate 2D window
%
% SYNOPSIS ImageProfile(imgH)
%
% INPUT	imgH	:	an image handle OR a figure handle
%
% c: 16/6/99	dT


% Check whether we got an image or a figure
if (strcmp(get(imgH,'Type'),'image'))
   imAxesH = get(imgH,'Parent');
   imFig = get(imAxesH,'Parent');
   img = get(imgH,'CData');
   figure(imFig);
else
   imFig=imgH;
   figure(imFig);
  	imAxesH = get(imFig,'CurrentAxes');
   imgH = findobj(imAxesH,'Type','image','Parent',imAxesH);
   img = get(imgH,'CData');
end;
% Dont do anything, if there's no image
if isempty(img)
   return;
end;

color = ['r' 'g' 'b' 'c' 'm' 'y'];		%Colortable for the line segments
imSize = size(img);
[xPts,yPts]=getline;
nPts= length(xPts);
profDat = [];
nSteps(1) =0;
%	Get color profile from line segments
for i = 1:(nPts-1)
   lineH(i)=line([xPts(i) ; xPts(i+1)],[yPts(i) ; yPts(i+1)],'Color',color(mod(i-1,6)+1));
   [xDum, YDum, proDum] = improfile(img,[xPts(i) ; xPts(i+1)],[yPts(i) ; yPts(i+1)]);
   profDat=[profDat;proDum];
   nSteps(i+1)=nSteps(i)+length(proDum)-1;
end;
% create new figure
proFig = figure;
set(proFig,'NumberTitle','off','Name','Color Profile');
set(proFig,'UserData',lineH);
set(proFig,'DeleteFcn','uiProfDeleteFcn')
% draw all plots in the same graph with different colors
hold on;	
for i = 1:(nPts-1)
   plot(nSteps(i):nSteps(i+1),profDat((nSteps(i)+1):(nSteps(i+1)+1)),'Color',color(mod(i-1,6)+1));
end;
hold off;
grid on;