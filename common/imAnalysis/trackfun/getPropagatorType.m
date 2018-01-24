function res = getPropagatorType()
%GETPROPAGATORTYPE query current propagator type
%
% SYNOPSIS res = getPropagatorType()
%
% INPUT : none
%
% OUTPUT res : string defining the current propagator

global prop__;

res = prop__.type;
