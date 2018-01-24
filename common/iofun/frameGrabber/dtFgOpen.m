function dtFgOpen
%DTFGOPEN opens the DataTranslation framegrabber for communications
%
% SYNOPSIS  dtFgOpen
% INPUT     none
% OUTPUT    none
%
% NOTE      Function refers to the global variable dtFgIsOpen__ 
%           which is set to 0 (false) when starting MATLAB through
%           $MATLAB/toolbox/local/startup.m

global dtFgIsOpen__;

if dtFgIsOpen__ > 0 
   return;
end;	
% call the mexFunction which opens the device
[status, dtFgIsOpen__ ] = mexDTFgHandler('open');
return;  