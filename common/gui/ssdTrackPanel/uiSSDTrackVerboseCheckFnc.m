function uiSSDTrackVerboseCheckFnc
% create function for the verbose checkbox in the SSD Track panel
global UISSDTRACKVERBOSE__;

pI = get(gcbf,'UserData');

val = get(gcbo,'Value');
pI.trackCmd.verbose = val;
UISSDTRACKVERBOSE__ = val;
set(gcbf,'UserData',pI);
