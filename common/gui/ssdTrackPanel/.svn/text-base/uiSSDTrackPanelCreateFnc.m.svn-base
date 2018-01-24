function uiSSDTrackPanelCreateFnc;
% create function for SSD Track panel
global ssdDfltsCmd__;
global ssdPlotLogFileName__;
global UISSDTRACKVERBOSE__;

set(gcbf,'Colormap',[]);
set(gcbf,'MinColormap',0);

UISSDTRACKVERBOSE__ = 0;

% initialize the panel info
pI.connectedView = [];
pI.stack = [];
pI.stackName = [];
pI.stackIndx = [];
pI.stackTot = [];
pI.nextIndx = [];
pI.nImgGrabbed = 0;
pI.tImg = [];
pI.tImgMaskType = 0;
pI.sImg = [];
pI.tPos = [];
pI.tPosOff = [];
pI.tDispOri = [];
pI.tPosRef = [0,0];
pI.pos0 = [];
pI.pos = [];
pI.posOff = [];
pI.patchDim = [];
pI.trackCmd = ssdDfltsCmd__;

set(gcbf,'UserData',pI);

% create a propagator 
initPropagator('ID');

% initialize the framegrabber board
openFg;

% delete the old plotting log file
if(exist(ssdPlotLogFileName__,'file'))
   delete(ssdPlotLogFileName__);
end;
