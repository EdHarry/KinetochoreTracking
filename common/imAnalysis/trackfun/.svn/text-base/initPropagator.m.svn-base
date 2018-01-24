function initPropagator(newType)
%INITPROPAGATOR initializes a coordinate propagator
%
% SYNOPSIS initPropagator(newType)
%
% INPUT newType : propagator type; the following types are implemented
%         ID -> position and shape left unchanged
%         IDNOSHAPE -> position unchanged, shape is returned as [1,0; 0,1]]
%         2PTOSCI -> oscillation between 2 points
%
% OUTPUT : none
%
% NOTE : setup of a global variable prop__

global prop__;

% check validity of type and set it 
switch newType,
case 'ID', prop__.type = newType;
case 'IDNOSHAPE', prop__.type = newType;
case '2PTOSCI', prop__.type = newType;
otherwise, prop__.type = 'NONE';
end;

prop__.currentPos = [];
prop__.currentShape = [];
prop__.nextPos = [];
prop__.nextShape = [];

prop__.params = [];
prop__.history = [];
