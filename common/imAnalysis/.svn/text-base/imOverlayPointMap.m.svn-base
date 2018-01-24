function imOverlayPointMap(map,fmt)
%IMOVERLAYPOINTMAP overlays an point (e.g. locmax) map to an image in a figure
%
% SYNOPSIS imOverlayPointMap(map,fmt)
%
% INPUT map : binary map with 1 at point location 0 otherwise 
%       fmt : plot format (see plot() for further information)
%             recommended either 'r.' or 'r+'
%
% SEE ALSO plot

indx = find(map);
cols = floor(indx/size(map,1))+1;
rows = rem(indx,size(map,1));

% an image is expected in the current figure
hold on;
plot(cols,rows,fmt);
hold off;
