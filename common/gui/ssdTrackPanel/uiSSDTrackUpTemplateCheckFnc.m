function uiSSDTrackUpTimgCheckFnc(type)
%UISSDTRACKUPTEMPLATECHECKFNC callback handling actions of 'update template' check box
%
% SYNOPSIS uiSSDTrackUpTimgCheckCreateFnc(create)
%
% INPUT type : 'create' -> callback during creation of the object
%              'pressed' -> callback after clicking

switch type,
case 'create', set(gcbo,'Value',0);
case 'pressed', setConsistentProp(get(gcbo,'Value'));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions 
function setConsistentProp(val)
% the update template function is not compatible 
% with the propagator setting 'ID'
if (val & isPropagatorType('ID'))
   changePropagatorType('IDNOSHAPE');
end;

