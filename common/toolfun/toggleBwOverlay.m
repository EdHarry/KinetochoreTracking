function toggleBwOverlay()
%TOGGLEBWOVERLAY toggles the B & W overlay status
%
% SYNOPSIS toggleBwOverlay
%
% INPUT none;
%
% OUTPUT none;

global bwOverlays__;

bwOverlays__ = mod(bwOverlays__ + 1,2);