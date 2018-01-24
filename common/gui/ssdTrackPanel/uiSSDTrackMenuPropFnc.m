function uiSSDTrackMenuPropFnc
%UISSDRACKMENUPROPFNC handles all the callbacks of the propagator menu

whoCalls = get(gcbo,'Tag');

switch whoCalls
case 'UISSDTRACKMENU_PROP', masterMenu;
case 'UISSDTRACKMENU_PROP_PROPS', propsMenu;
case 'UISSDTRACKMENU_PROP_ID', newProp('ID');
case 'UISSDTRACKMENU_PROP_IDNOSHAPE', newProp('IDNOSHAPE'); 
case 'UISSDTRACKMENU_PROP_2PTOSCI', newProp('2PTOSCI');     
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  local functions

function masterMenu
subMenus = get(gcbo,'Children');

for(i=1:length(subMenus))
   set(subMenus(i),'Checked','off');
end;

set(findobj(gcbo,'Tag','UISSDTRACKMENU_PROP_PROPS'),'Enable','off');

if(isPropagatorType('ID'))
   set(findobj(gcbo,'Tag','UISSDTRACKMENU_PROP_ID'),'Checked','on');
   set(findobj(gcbo,'Tag','UISSDTRACKMENU_PROP_PROPS'),'Enable','off');
   return;
end;
if(isPropagatorType('IDNOSHAPE'))
   set(findobj(gcbo,'Tag','UISSDTRACKMENU_PROP_IDNOSHAPE'),'Checked','on');
   set(findobj(gcbo,'Tag','UISSDTRACKMENU_PROP_PROPS'),'Enable','off');
   return;
end;
if(isPropagatorType('2PTOSCI'))
   set(findobj(gcbo,'Tag','UISSDTRACKMENU_PROP_2PTOSCI'),'Checked','on');
   set(findobj(gcbo,'Tag','UISSDTRACKMENU_PROP_PROPS'),'Enable','on');
   return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function propsMenu

if(isPropagatorType('2PTOSCI'))
   uiSSDTrackPanelSet2PTOSCIparams;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newProp(type)

% with the currently available propagators 'ID' is the one 
% which also propagates the shape matrix. This, however, is not
% possible as long as the checkbox 'update template' is activated
if(strcmp(type,'ID') & get(findobj(gcbf,'Tag','UISSDTRACKUPDATETIMGCHECK'),'Value'))
   return;
else
   set(gcbo,'Checked','on');
   changePropagatorType(type);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



