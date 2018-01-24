function [pos, shape] = getPropagatorNextPos()
%GETPROPAGATORNEXTPOS returns the propagated next position.
%
% SYNOPSIS : [pos, shape] = getPropagatorNextPos()
%
% INPUT : none 
%
% OUTPUT : pos propagated position in a left handed coordinate system with origin (1,1)

global prop__;

if(isempty(prop__.currentPos))
   return;
end;

switch prop__.type
case 'ID', setIDPropNextPos;
case 'IDNOSHAPE', setIDPropNextPos;
case '2PTOSCI', set2PTOSCIPropNextPos;
end;

pos = prop__.nextPos;
shape = prop__.nextShape;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  local functions

function setIDPropNextPos()
% service function for the ID propagator

global prop__;

% just set the next position as the current
prop__.nextPos = prop__.currentPos;

switch prop__.type
case 'ID', prop__.nextShape = prop__.currentShape;
case 'IDNOSHAPE', prop__.nextShape = [];
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set2PTOSCIPropNextPos()
% service function which sets the next position for the coordinate propagator

global prop__;

if(isempty(prop__.params))
   % no parameters have been set so far
   prop__.nextPos = prop__.currentPos;
else
   newPos(1) = cos(prop__.params(1))* prop__.params(2);
   newPos(2) = sin(prop__.params(1))* prop__.params(2);
   prop__.nextPos = prop__.currentPos + newPos;
end;
