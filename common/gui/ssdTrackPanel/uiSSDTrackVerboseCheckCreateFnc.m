function uiSSDTrackVerboseCheckCreateFnc
% create function for the verbose checkbox in the SSD Track panel

global UISSDTRACKVERBOSE__;

pI = get(gcbf,'UserData');

pI.trackCmd.verbose = UISSDTRACKVERBOSE__;
set(gcbf,'UserData',pI);

set(gcbo,'Value',UISSDTRACKVERBOSE__);
