function uiSSDTrackMenuDelStackFnc
% call back for the Delete Stack Menu of the SSD Track panel
pI = get(gcbf,'UserData');

pI.stack = [];
pI.stackName = [];
pI.stackIndx = [];
pI.stackTot = [];
pI.nextIndx = [];
pI.nImgGrabbed = 0;
pI.tImg = [];
pI.sImg = [];
pI.tPos = [];
pI.tPosOff = [];
pI.tDispOri = [];
pI.pos0 = [];
pI.pos = [];
pI.posOff = [];
pI.patchDim = [];

set(findobj(gcbf,'Tag','UISSDTRACKNEXTBUTTON'),'Enable','off');
set(findobj(gcbf,'Tag','UISSDTRACKGOBUTTON'),'Enable','off');
set(findobj(gcbf,'Tag','UISSDTRACKWRAPCHECK'),'Enable','off');
set(findobj(gcbf,'Tag','UISSDTRACKWRAPTEXT'),'Enable','off');
set(findobj(gcbf,'Tag','UISSDTRACKFROMLOGCHECK'),'Enable','on');
if(~isReadyFg )
   set(findobj(gcbf,'Tag','UISSDTRACKINITBUTTON'),'Enable','off');
end;



set(gcbf,'UserData',pI);

% make text of the view panel invisible
pI.connectedView = queryUiViewPanel(pI.connectedView);
if(~isempty(pI.connectedView))
   textH = findobj(pI.connectedView,'Type','uicontrol','Tag','UIVIEWTEXT');
   set(textH,'Visible','off');
end;

