function dbstopIfLasterror
%DBSTOPIFLASTERROR sets a breakpoint where the last error occured
%
% SYNOPSIS dbstopIfLasterror
%
% REMARKS  the function will query lasterror and stop the next time it
%          finds this particular error.
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: jdorn
% DATE: 21-Feb-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% query lasterror
err = lasterror;

% set breakpoint
if ~isempty(err)
    eval(sprintf('dbstop if caught error %s', err.identifier));
    % talk to user
    disp(sprintf('Matlab will break at error %s',err.message));
end
