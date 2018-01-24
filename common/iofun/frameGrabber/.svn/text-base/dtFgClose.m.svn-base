function dtFgClose
% Closes communications to the DataTranslation framegrabber
%
% SYNOPSIS  DTFgClose
% INPUT     none
% OUTPUT    none
%
% NOTE      Function refers to the global variable dtFgIsOpen__ 
%           which is set to 0 (false) when starting MATLAB through 
%           $MATLAB/toolbox/local/startup.m

global dtFgIsOpen__;

if dtFgIsOpen__ == 0 
   return;
end;	
% call the mexFunction which closes the device
[status, dtFgIsOpen__ ] = mexDTFgHandler('close');
return;  