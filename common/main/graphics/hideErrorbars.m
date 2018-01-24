function hideErrorbars(figH)
%HIDEERRORBARS hides all errorbars in a figure
%
% SYNOPSIS hideErrorbars(figH)
%
% INPUT    figH handle of the figure you want to turn the errorbars off
%           complement to unhideErrorbars). Errorbars have to have been
%           drawn with myErrorbar.m or they have to be tagged 'errorBar'
%
% c: 05/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% is there a figure
if nargin == 0 | isempty(figH)
    figH = get(0,'CurrentFigure');
    if isempty(figH)
        error('no figure open')
    end
end
    

% find lineHandles in figure
lineHandles = findall(figH,'Type','line');

for lh = lineHandles'
    if findstr(get(lh,'Tag'),'errorBar')
        set(lh,'Visible','off')
    end
end