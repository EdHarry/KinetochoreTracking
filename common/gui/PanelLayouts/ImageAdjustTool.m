function fig = ScaleOptions1()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.

load ScaleOptions1

h0 = figure('Color',[0.8 0.8 0.8], ...
	'Colormap',mat0, ...
	'MenuBar','none', ...
	'Name','Gamma Correction', ...
	'NumberTitle','off', ...
	'PointerShapeCData',mat1, ...
	'Position',[513 328 285 360], ...
	'Resize','off', ...
	'Tag','IMAGE_ADJUST_TOOL', ...
	'UserData',mat2);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','ImAdjustToolCB(get(gcf,''Userdata''))', ...
	'ListboxTop',0, ...
	'Position',[157.5 11.25 45 15], ...
	'String','Close', ...
	'Tag','CLOSE');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[123.75 236.25 33.75 15], ...
	'String','0', ...
	'Style','edit', ...
	'Tag','IN1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[168.75 236.25 33.75 15], ...
	'String','1', ...
	'Style','edit', ...
	'Tag','IN2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[11.25 236.25 90 11.25], ...
	'String','Input intensity interval:', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[11.25 213.75 90 11.25], ...
	'String','Output intensity interval:', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[168.75 213.75 33.75 15], ...
	'String','1', ...
	'Style','edit', ...
	'Tag','OUT2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[123.75 213.75 33.75 15], ...
	'String','0', ...
	'Style','edit', ...
	'Tag','OUT1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','ImAdjustToolCB(get(gcf,''Userdata''))', ...
	'ListboxTop',0, ...
	'Position',[50 191.25 45 15], ...
	'String','0.26874', ...
	'Style','edit', ...
	'Tag','GAMMAEDIT');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'Callback','ImAdjustToolCB(get(gcf,''Userdata''))', ...
	'ListboxTop',0, ...
	'Max',2, ...
	'Position',[101.25 191.25 101.25 15], ...
	'Style','slider', ...
	'Tag','GAMMA', ...
	'Value',0.7200000000000002);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[11.25 191.25 30 11.25], ...
	'String','Gamma:', ...
	'Style','text', ...
	'Tag','StaticText2');
h1 = axes('Parent',h0, ...
	'Units','pixels', ...
	'Box','on', ...
	'CameraUpVector',[0 1 0], ...
	'Color',[1 1 1], ...
	'ColorOrder',mat3, ...
	'Position',[31 61 240 180], ...
	'Tag','Axes1', ...
	'XColor',[0 0 0], ...
	'XGrid','on', ...
	'XLimMode','manual', ...
	'YColor',[0 0 0], ...
	'YGrid','on', ...
	'YLimMode','manual', ...
	'ZColor',[0 0 0]);
h2 = line('Parent',h1, ...
	'Color',[0 0 1], ...
	'Tag','Axes1Line3', ...
	'XData',mat4, ...
	'YData',mat5);
h2 = line('Parent',h1, ...
	'Color',[0 0 1], ...
	'Tag','Axes1Line2', ...
	'XData',[0 0], ...
	'YData',[0 0]);
h2 = line('Parent',h1, ...
	'Color',[0 0 1], ...
	'Tag','Axes1Line1', ...
	'XData',[1 1], ...
	'YData',[1 1]);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[0.497907949790795 -0.1340782122905029 17.32050807568877], ...
	'Tag','Axes1Text4', ...
	'VerticalAlignment','cap');
set(get(h2,'Parent'),'XLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[-0.1213389121338912 0.4916201117318435 17.32050807568877], ...
	'Rotation',90, ...
	'Tag','Axes1Text3', ...
	'VerticalAlignment','baseline');
set(get(h2,'Parent'),'YLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','right', ...
	'Position',[-0.1297071129707113 1.664804469273743 17.32050807568877], ...
	'Tag','Axes1Text2', ...
	'Visible','off');
set(get(h2,'Parent'),'ZLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',mat6, ...
	'Tag','Axes1Text1', ...
	'VerticalAlignment','bottom');
set(get(h2,'Parent'),'Title',h2);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','ImAdjustToolCB(get(gcf,''Userdata''))', ...
	'ListboxTop',0, ...
	'Position',[22.5 11.25 45 15], ...
	'String','Reset', ...
	'Tag','RESET');
if nargout > 0, fig = h0; end
