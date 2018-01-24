function stk2tif(filename)

% stk2tif converts STK to TIF using metaTiffRead
%
% SYNOPSIS   stk2tif(filename)
%
% INPUT      filename : string containing the common name of the sequence of .tif files
%
% OUTPUT     none     :    
% 
% REMARKS    files are written to disk into user-selected directory
%       
%
%
% DEPENDENCES   stk2tif uses {metaTiffRead}
%               stk2tif is used by {}
%
% Alexandre Matov, November 7th, 2002

if nargin~=1
    error('Please enter a valid (common) file name for the output files');
end
if isempty(filename)
    error('Please enter a valid (common) file name for the output files');
end

try
    [S,n]=metaTiffRead;
catch
    disp('Interrupted by user.');
    return
end
L=length(num2str(n)); 
strg=sprintf('%%.%dd',L); % Creates the format string for the numerical indexes

% Select a directory where the output .tif files will be written
path=uigetdir('','Select output directory');
if path==0
    disp('Aborting...');
    return
end    

h=waitbar(0,'Please wait! Writing .tif files');
for i=1:n
    indxStr=sprintf(strg,i);
    imwrite(S(i).data,[path,filesep,filename,indxStr,'.tif']);
    waitbar(i/n,h);
end
close(h);