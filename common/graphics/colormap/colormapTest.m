function fh = colormapTest(cmap,name)
%COLORMAPTEST shows test plots for colormaps
%
% SYNOPSIS: fh = colormapTest(cmap,name)
%
% INPUT cmap: n-by-3 RGB colormap
%       name: name of colormap
%
% OUTPUT fh: handle to test figure
%
% REMARKS
%
% created with MATLAB ver.: 7.3.0.267 (R2006b) on Windows_NT
%
% created by: jonas
% DATE: 02-Apr-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check colormap
if nargin == 0 
    error('please supply a colormap!')
end
% check size, orientation
cmap = returnRightVector(cmap, 3, 'c');

% make test images
% 1) 2 Gaussians moving together
% 2) peaks function
gaussX = repmat(-10:0.2:10,101,1);
gaussY = repmat((0:100)',1,101);
gaussZ = exp(-(gaussX-6+0.06*gaussY).^2./2) + ...
    exp(-(gaussX+6-0.06*gaussY).^2./2) + ...
    0.5*exp(-(gaussY-40+gaussX).^2./200);

peaksZ = peaks(100);


% plot. Top: Gauss, Bottom: Peaks
if nargin < 2 || isempty(name)
    handle = figure('Name',sprintf('Test of colormap %s',inputname(1)));
else
    handle = figure('Name',sprintf('Test of colormap %s',name));
end


% plot grey with light - can't do grey!!
ah = subplot(2,3,1);
surf(gaussX,gaussY,gaussZ,'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
set(ah,'xGrid','off','yGrid','off','zGrid','off','Color','none','xTick',[],'yTick',[],'zTick',[])
%colormap gray
light('Position',[-10,100,2])
light('Position',[-10,0,2])
material shiny

% plot 3d with colormap
ah = subplot(2,3,2);
surf(gaussX,gaussY,gaussZ,'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
set(ah,'xGrid','off','yGrid','off','zGrid','off','Color','none','xTick',[],'yTick',[],'zTick',[])
colormap(cmap)

% plot 2d - with colorbar
ah = subplot(2,3,3);
contourf(gaussX,gaussY,gaussZ,'LineStyle','none','LevelList',linspace(nanmin(gaussZ(:)),...
    nanmax(gaussZ(:)),100));
set(ah,'xGrid','off','yGrid','off','zGrid','off','Color','none','xTick',[],'yTick',[],'zTick',[])
colormap(cmap)
axis square
colorbar('peer',ah)

% do the same for peaks
ah = subplot(2,3,4);
surf(peaksZ,'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
set(ah,'xGrid','off','yGrid','off','zGrid','off','Color','none',...
    'xTick',[],'yTick',[],'zTick',[],'clim',[-9,9])
%colormap gray
light('Position',[-100,100,10])
light('Position',[-100,0,10])
material shiny

% plot 3d with colormap
ah = subplot(2,3,5);
surf(peaksZ,'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
set(ah,'xGrid','off','yGrid','off','zGrid','off','Color','none',...
    'xTick',[],'yTick',[],'zTick',[],'clim',[-9,9])
colormap(cmap)

% plot 2d - with colorbar
ah = subplot(2,3,6);
contourf(peaksZ,'LineStyle','none','LevelList',linspace(nanmin(peaksZ(:)),...
    nanmax(peaksZ(:)),100));
set(ah,'xGrid','off','yGrid','off','zGrid','off','Color','none',...
    'xTick',[],'yTick',[],'zTick',[],'clim',[-9,9])
colormap(cmap)
axis square
colorbar('peer',ah)

% return figure handle if requested
if nargout > 0
    fh = handle;
end
