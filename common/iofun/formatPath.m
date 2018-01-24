function [newPath]=formatPath(oldPath)
% FORMATPATH converts between linux and window paths and vice versa
%
% this function attempts to create a directory name from the input path as
% well as the currently working directory. if the path generated is not an
% existing directory, the function prompts the user to select another
% directory. the user should then select any directory above or below the
% input directory, or the input directory itself (making sure it does
% in fact exist).  the reason this may occur is the current working
% directory may not be pointing to the same server where the input
% directory exists.
%
% Kathryn Applegate 2008


% switch direction of fileseps
if ispc
    temp=strrep(oldPath, '/', '\');       
else
    temp=strrep(oldPath, '\', '/');
end
if isequal(temp,oldPath)
    % OS didn't change, nothing to do.
    newPath=oldPath;
    return
end

% check to make sure the input path doesn't contain white space
whiteSpaceIdx=regexp(temp,['\s'],'start')';
if ~isempty(whiteSpaceIdx)
    error('formatPath: input directory name must not include spaces')
end

doneFlag=0;
% look at current directory
currentDir=pwd;
% find oldPath's filesep locations
tempFilesepIdx=strfind(temp,filesep);

while doneFlag==0
    % find current directory's filesep locations
    currentDirFilesepIdx=strfind(currentDir,filesep);
    
    if ispc
        % replace, for example, '/mnt/fsm' with 'S:'
        newPath=[currentDir(1:currentDirFilesepIdx(1)-1) temp(tempFilesepIdx(3):end)];
    else
        % OS is linux
        % there should be at least 3 fileseps for a linux directory
        % (i.e. '/mnt/fsm/'), but there may only be 1 or 2.  in this case
        % the name will not yield a real directory, but the user will be
        % directed to select one on the path. here we take care of indexing
        % problem in case of 1 or 2.
        num=min(length(currentDirFilesepIdx),3);
        % replace, for example, 'S:\' with '/mnt/fsm'
        newPath=[currentDir(1:currentDirFilesepIdx(num)-1) temp(tempFilesepIdx(1):end)];
    end

    % check if the created path is actually a directory. if not,the root
    % was wrong. ask the user to select a new one.
    if isdir(newPath)
        doneFlag=1;
    else
        currentDir=uigetdir(currentDir,'Select any directory above input directory');
    end
end

