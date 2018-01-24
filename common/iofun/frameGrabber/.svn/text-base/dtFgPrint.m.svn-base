function dtFgPrint
% Prints the settings of DataTranslation framegrabber
%
% SYNOPSIS  dtFgPrint
% INPUT     none
% OUTPUT    none
%
% NOTE      Function refers to the global variable dtFgIsOpen__
%           and runs only if this is set to TRUE

global dtFgIsOpen__;

if dtFgIsOpen__ == 0 
   return;
end;	
% call the mexFunction which prints the device info
[status, dtFgIsOpen__] = mexDTFgHandler('print');
return;  