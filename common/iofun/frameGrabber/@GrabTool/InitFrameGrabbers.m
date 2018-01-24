function [fgTypeList,nrOfGrabbers] = InitFrameGrabbers(gbT)
%GRABTOOL/INITFRAMEGRABBERS checks and initalizes installed framegrabbers
%
% SYNOPSIS OpenGrabPanel(gbT)
%
% INPUT  gbT: an grabTool object
%
% OUTPUT fgTypeList: List of connected framegrabber types
%        nrOfGrabbers: vector of the # of installed grabbers
%                      in the same order as fgTypeList
% c: 24/8/99	dT

fgTypeList={'None' '2 x PP Mono'};
nrOfGrabbers=2;