function ans = dtFgIsOpen
%DTFGISOPEN Checks the status of the DataTranslation framegrabber
%
% SYNOPSIS  dtFgIsOpen
% INPUT     none
% OUTPUT    ans : >0 if open
%                 =0 if not open
%
% NOTE      Function refers to the global variable dtFgIsOpen__ 
%           which is set to 0 (false) when starting MATLAB through
%           $MATLAB/toolbox/local/startup.m

global dtFgIsOpen__;
ans = (dtFgIsOpen__ > 0);