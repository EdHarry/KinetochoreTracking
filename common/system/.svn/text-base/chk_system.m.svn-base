function ans=chk_system()
% CHK_SYSTEM get the operating system
% 
%             
% SYNOPSIS        ans=chk_system()
%
% INPUT           
% 
% OUTPUT          ans: 1: is Linux
%                      2: is Unix
%                      3: should be Windows
%                           
% DEPENDENCES     createLogFile { searchFiles
%                                 nowstring  
%                               }
% Matthias Machacek 29/10/03

[stat, sys]=system('uname');
if  strncmp(sys,'Linux',3)
   ans=1;       
elseif  strncmp(sys,'Unix',3)
   ans=2;     
else
   %this must be probably a windows system
   ans=3;
end