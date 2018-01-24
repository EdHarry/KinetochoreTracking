%DBCC clears classes and sets dbstop if error
%
% SYNOPSIS: dbcc
%
% INPUT 
%
% OUTPUT 
%
% REMARKS
%
% created with MATLAB ver.: 7.7.0.2162 (R2008b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 01-Oct-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dbs = dbstatus;
save('dbstatus.mat','dbs')
clear classes
load('dbstatus.mat')
for i=1:length(dbs)
    try
        dbstop(dbs(i));
    catch
        fprintf('could not restore dbstop %s\n',dbs(i).name);
    end
end
delete('dbstatus.mat')
clear dbs i