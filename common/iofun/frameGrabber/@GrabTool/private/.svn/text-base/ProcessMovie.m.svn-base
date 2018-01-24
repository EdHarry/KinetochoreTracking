function gbT = ProcessMovie(gbT)
%GRABTOOL/PROCESSMOVIE
% This method creates/shows a movie depending on wether
% the stack has already been converted or not

if(isempty(gbT.movie.fps) | (gbT.movie.fps<= 0))
   gbT.movie.fps = 1;   
end;


% check whether there is a view 
if( isempty(gbT.viewPanelH) )
   gbT.viewPanelH = uiViewPanel;
else
   gbT.viewPanelH = ...
      queryUiViewPanel(gbT.viewPanelH);
   if( isempty(gbT.viewPanelH) )
      gbT.viewPanelH = uiViewPanel;
   end;
end;

% get the axes handle in this view
axesH = get(gbT.viewPanelH,'CurrentAxes');

% check whether the movie has already been generated
if( isempty(gbT.movie.data))
   [gbT.movie.data, gbT.movie.cMap] = ...
      stack2movie(gbT.data,axesH,128,1,gbT.movie.fps);
else
   figure(gbT.viewPanelH);
   colormap(gbT.movie.cMap);
   if(isempty(gbT.movie.repeat) | (gbT.movie.repeat <= 0))
      gbT.movie.repeat=1;
   end;
   movie(gbT.movie.data,gbT.movie.repeat,gbT.movie.fps);
end;
