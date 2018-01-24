function [listOfFunctions, toolboxes, neededFunctions] = requiredFunctions(functionList)
%REQUIREDFUNCTIONS returns functions required to run certain code
%
% SYNOPSIS: [listOfFunctions, toolboxes] = requiredFunctions(functionList)
%
% INPUT functionList : cell array containing functions to be tested (or
%               string with single function). 
%               Suggestion: search for all the functions in a directory
%               with searchFiles('.m$','','ask'), and pass the first column
%               of the output to requiredFunctions.
%
% OUTPUT listOfFunctions : list of functions from MatlabRoot needed to run
%               the specified functions. Does not contain operator
%               functions (functions that are specifically defined for
%               certain variable classes other than double).
%               Note that the input functions will also show up here.
%		 toolboxes : names of required Matlab toolboxes 
%        neededFunctions : result of depfun(functionList{:})
%
%
% REMARKS - Running this function can take a while because it calls depfun.
%         - Use the output with copyfile to copy functions to an output
%           directory, if your goal is to share code.
%
% created with MATLAB ver.: 7.5.0.342 (R2007b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 08-Jan-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% modify input
if ~iscell(functionList)
    functionList = {functionList};
end


% check for list of functions
neededFunctions = depfun(functionList{:},'-quiet');

% find matlab's own functions
matlabIdx = strmatch(matlabroot,neededFunctions);

% get toolboxes
[dummy,dummy,dummy,dummy,tokens] = regexp(neededFunctions(matlabIdx),'\\toolbox\\(\w+)\\');
aa = cat(1,tokens{:});aa = cat(1,aa{:});
aaa = strvcat(aa{:}); %#ok<VCAT>
toolboxes = unique(aaa,'rows');

% create list of own functions 
listOfFunctions = neededFunctions;
listOfFunctions(matlabIdx) = [];



% remove operators from list
opIdx = regexp(listOfFunctions,'\\@.?int');
opIdxL = ~cellfun(@isempty,opIdx);
listOfFunctions(opIdxL) = [];

% for listOfFunctions: Check whether there are gui-figs that are required
for f = 1:length(listOfFunctions)
    [pname,fname] = fileparts(listOfFunctions{f});
    guiName = fullfile(pname,[fname,'.fig']);
    if exist(guiName,'file')
        listOfFunctions{end+1} = guiName; %#ok<AGROW>
    end
end
%... and sort
listOfFunctions = sort(listOfFunctions);

