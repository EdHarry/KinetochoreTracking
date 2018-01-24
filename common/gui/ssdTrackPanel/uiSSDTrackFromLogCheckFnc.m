function uiSSDTrackFromLogCheckFnc(type)
%UISSDTRACKFROMLOGCHECKFNC callback handling actions of 'update template' check box
%
% SYNOPSIS uiSSDTrackFromLogCheckFnc(type)
%
% INPUT type : 'create' -> callback during creation of the object
%

switch type,
case 'create', set(gcbo,'Value',0); set(gcbo,'Enable','off');
end;

