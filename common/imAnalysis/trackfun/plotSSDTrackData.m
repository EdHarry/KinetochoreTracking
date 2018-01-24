function newFigH = plotSSDTrackData(data,refIndx,opt,figH)
%PLOTSSDTRACKDATA plots the results in a log file from SSD tracking
%
% SYNOPSIS newFigH = plotSSDTrackData(data,refIdx,opt,figH)
%
% INPUT data : data matrix found in the SSD Tratck log files
%       refIdx : index of a reference point (0 means plotting absolute coordinates)
%       opt : structure with plot options (all the fields are optional)
%             *.hold : 1 if an old plot in figH should be hold
%             *.linetype  : linetype defition according to the convention in 
%                           MATLAB plot()
%             *.axeslim : limit of the axes [x1min,x1max,x2min,x2max]; if only
%                         one pair is given both axis get the same limits
%             *.startmark : special mark for the first data point (see MATLAB plot() ) 
%             *.endmark : special mark for the last data point (see MATLAB plot() )
%             *.pixsze : pixel size in um [s1, s2]; if this field is set the plot 
%                        is converted to microns; if only one value is given, 
%                        square pixels are assumed.
%       figH : (optional) figure handle 
%
% OUTPUT newFigH : new figure handle
%

holdOpt = 0;

if(isfield(opt,'linetype'))
   plotopt = opt.linetype;
else
   plotopt = 'r-';
end;

x1Text = 'x_1 position [pixel]';
x2Text = 'x_2 position [pixel]';

% check validity of the input figure handle
if(nargin == 4)
   if(isempty(findobj(figH,'Type','figure')))
      newFigH = figure;
   else
      newFigH = figure(figH);
      if(isfield(opt,'hold'))
         holdOpt = 1;
      end;
   end;
else
   newFigH = figure;
end;

% distinguish between different file formats 
switch size(data,2)
case 7, pos = [data(:,2),data(:,3)];
case 11, pos = [data(:,2),data(:,3)];
case 8, pos = [data(:,3),data(:,4)];
case 12, pos = [data(:,3),data(:,4)];
otherwise error('unsupported file format');
end;

% eventually scale the data
if(isfield(opt,'pixsze'))
   fact(1) = opt.pixsze(1);
   if(length(opt.pixsze)==2)
      fact(2) = opt.pixsze(2);
   else
      fact(2) = opt.pixsze(1);
   end;
   pos(:,1) = pos(:,1).*fact(1);
   pos(:,2) = pos(:,2).*fact(2);
   x1Text = 'x_1 position [um]';
   x2Text = 'x_2 position [um]';
end;

if((refIndx > 0)& (refIndx <= size(data,1)))
   pos(:,1) = pos(:,1) - pos(refIndx,1);
   pos(:,2) = pos(:,2) - pos(refIndx,2);
end;

% plot the data
if(holdOpt) 
   hold on;
end;
lH = plot(pos(:,1),pos(:,2),plotopt);
set(lH,'LineWidth',1.5);
if(holdOpt) 
   hold off;
end;

% set start mark
if(isfield(opt,'startmark'))
   hold on;
   lH = plot(pos(1,1),pos(1,2),opt.startmark);
   set(lH,'MarkerSize',25);
   hold off;
end;

% set end mark'
endIndx = size(pos,1);
if(isfield(opt,'endmark'))
   hold on;
   lH = plot(pos(endIndx,1),pos(endIndx,2),opt.endmark);
   set(lH,'MarkerSize',25);
   hold off;
end;

% make aspect axes aspect ratio equal
set(gca,'PlotBoxAspectRatioMode','manual');
set(gca,'PlotBoxAspectRatio',[1,1,1]);

% set the axis limits
if(isfield(opt,'axeslim'))
   set(gca,'XLim',opt.axeslim(1:2));
   if(length(opt.axeslim) == 4)
      set(gca,'YLim',opt.axeslim(3:4));
   else
      set(gca,'YLim',opt.axeslim(1:2));
   end;
end;

% set the axes labels
xlabel(x1Text);
ylabel(x2Text);

  
   

      

