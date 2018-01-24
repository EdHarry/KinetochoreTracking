function ret = getPropagatorParams
%GETPROPAGATORPARAMS gets the parameter vector of the coordinate propagator
%
% SYNOPSIS ret = getPropagatorParams();
% 
% INPUT none
%
% OUTPUT ret : parameter vector
%              dependent on the propagator type the vector components have different
%              meanings:
%       '2PTOSCI' : [azimuth [0,2pi], excursion [0,inf], nDataPoints]
%
%       if no parameter vector is defined the return value is an empty vector

global prop__;

if(isempty(prop__))
   ret = [];
   return;
end;

ret = prop__.params;
return;
   
