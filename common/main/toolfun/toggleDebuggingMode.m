function toggleDebuggingMode()
%TOGGLEDEBUGGINGMODE toggles the debuggingMode status for G.D. functions
%
% SYNOPSIS toggleDebuggingMode
%
% INPUT none;
%
% OUTPUT none;

global debuggingMode__;

debuggingMode__ = mod(debuggingMode__ + 1,2);