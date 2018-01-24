function setPropagatorPosEstimate(pos,shape)
%SETPROPAGATORPOSESTIMATE enters a new position measurement to the propogator.
% Dependent on the propagator type the current position will be updated
% accordingly
%
% SYNOPSIS setPropagatorPosEstimate(pos,shape)
%
% INPUT pos : position estimate in a left handed coordinate system with origin (1,1)
%       shape : [mx, sx; sy, my] (optional for those propagators which 
%                                 support shape propagation)
%
% OUTPUT : none


global prop__;

switch prop__.type
case 'ID', 
   if nargin == 1
      setIDPropPosEstimate(pos,[1,0;0,1]);
   else
      setIDPropPosEstimate(pos,shape);
   end;	
case 'IDNOSHAPE', setIDPropPosEstimate(pos,[1,0;0,1]);
case '2PTOSCI', set2PTOSCIPropPosEstimate(pos)
end;

   