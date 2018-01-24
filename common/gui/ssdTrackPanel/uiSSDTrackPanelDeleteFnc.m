function uiSSDTrackPanelDeleteFnc;
% callback for SSD track panel
global UISSDTRACKVERBOSE__;

% check if there is a an acquisition panel
acqH = findobj(0,'Type','figure','Tag','UIACQPANEL');

if(isReadyFg & isempty(acqH))
   % close the framegrabber connection properly
	closeFg;
end;


% delete the propagator variable 
deletePropagator;

% delete global variables associated with the GUI
clear global UISSDTRACKVERBOSE__;
