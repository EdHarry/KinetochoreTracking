function changePropagatorType(newType)
%CHANGEPROPAGATORTYPE changes the type of coordinate propagation
%
% SYNOPSIS changePropagatorType(newType)
%
% INPUT newType : new type of coordinate propagator
%                 the following types are implemented:
%                 ID -> position from instant to instant unchanged
%                 2PTOSCI -> oscillation between 2 points
%
% OUTPUT none

global prop__;

switch newType,
case 'ID', 
   prop__.type = newType;
   prop__.params = [];
case 'IDNOSHAPE',
   prop__.type = newType;
   prop__.params = [];
case '2PTOSCI', 
   prop__.type = newType;
   prop__.params = [];
otherwise, prop__.type = 'NONE';
   prop__.currentPos = []
   prop__.nextPos = [];
   prop__.params = [];
   prop__.history = [];
end;

   