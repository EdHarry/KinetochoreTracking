function res = isPropagatorType(type)
%ISPROPAGATORTYPE checks type of coordinate propagator
%
% SYNOPSIS res = isPropagatorType(type)
%
% INPUT type : string with type
%
% OUTPUT res : 1 if propagator is of type 'type'
%              0 otherwise

global prop__;

if(isempty(prop__))
   res = 0;
   return;
end;

if(strcmp(prop__.type,type))
   res = 1;
   return;
else
   res = 0;
end;
