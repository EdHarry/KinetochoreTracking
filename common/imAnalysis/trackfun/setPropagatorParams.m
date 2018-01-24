function setPropagatorParams(params)
%SETPROPAGATORPARAMS sets the parameter vector of the coordinate propagator.
% This function is only relevant for parametric propagators
%
% SYNOPSIS setPropagatorParams(params)
% 
% INPUT params : parameter vector; 
%                dependent on the propagator type the vector components have different
%                meanings:
%       '2PTOSCI' : [azimuth [0,2pi], excursion [0,inf], nDataPoints]
%
% OUTPUT none
%

global prop__;

switch(prop__.type)
case '2PTOSCI', 
   if(length(params) == 3)
      prop__.params = params;
   else 
      % invalid parameter vector 
      changePropagatorType('ID');
   end;   
otherwise,       
   changePropagatorType('ID');
end;
