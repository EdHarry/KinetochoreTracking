function uiViewPanelCreateFnc
% creation callback for the view panel

%set(gcbf,'MinColormap',128);

% cascades the new view Panel with respect to existing ones

handle = findobj('Type','figure','Tag','UIVIEWPANEL');
[nHandles,iDum] = size(handle);
if(nHandles > 1)
   pos = get(handle(2),'Position');
   pos(1) = pos(1) + 25;
   pos(2) = pos(2) - 25;
   set(gcbf,'Position',pos);
end;

