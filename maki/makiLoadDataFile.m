function dataStruct = makiLoadDataFile(serverType,dataFileName)
%MAKILOADDATAFILE loads a maki-data file from disk
%
% SYNOPSIS: dataStruct = makiLoadDataFile(serverType,dataFileName)
%
% INPUT     serverType: 'TEST', 'HERCULES', 'DANUSER', 'MERALDI', 'SWEDLOW'
%                       or 'MCAINSH'
%           dataFileName (opt): [pathName,filesep,filename] of dataFile. If
%           omitted, or if only a pathName is specified, the file can be
%           selected interactively
%
% OUTPUT dataStruct: data structure as described in makiMakeDataStruct
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 28-Jun-2007
%
% Edit by EHarry April 2012, now with version control on each individual
% field
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=================
%% CHECK INPUT
%=================

guiLoad = true;
guiPath = pwd;

if nargin < 2 || isempty(dataFileName)
    % guiLoad
    if nargin > 0 && ~isempty(serverType)
        guiPath = makiPathDef(['$' serverType],serverType);
    end
else
    % check for file, then for path, then exit with error
    if exist(dataFileName,'file')
        % we have a dataFile or a path
        [dataFilePath,part1,part2] = fileparts(dataFileName);
        if ~isempty(part2)
            % it's a data file (allow for stuff like .mat_old)
            dataFileName = [part1,part2];
            guiLoad = false;
        else
            % it's a path. fileParts will consider the last directory as
            % fileName
            guiPath = fullfile(dataFilePath,part1);
        end
    else
        error('file or path %s not found',dataFileName);
    end
    
end


%=================
%% GUILOAD
%=================

if guiLoad
    % there is a bug in Matlab 2007a that doesn't allow the use of filter
    % names
    %     oldDir = cd(guiPath);
    %      [dataFileName, dataFilePath] = uigetfile(...
    %          {'*-makiData-*','dataFiles';'*','all files'},'Please select makiDataFile');
    %      cd(oldDir);
    oldDir = cd(guiPath);
    [dataFileName, dataFilePath] = uigetfile(...
        '*-makiData-*','Please select makiDataFile');
    cd(oldDir);
    
    if dataFileName == 0
        error('no dataFile loaded')
    end
end

%================
%% LOAD DATA
%================

% load data file
load(fullfile(dataFilePath,dataFileName));
% write new path
dataStruct.dataFilePath = dataFilePath;
dataStruct.dataFileName = dataFileName;

% interpret rawMoviePath
if isempty(dataStruct.rawMoviePath)
    dataStruct.rawMoviePath = dataStruct.dataFilePath;
else
    % remove identifier
    dataStruct.rawMoviePath = makiPathDef(dataStruct.rawMoviePath,serverType);
end

% dependencies list count
depListC = 0;
namedFieldsC = 0;

% loop through dataStruct and read individual files
fn = fieldnames(dataStruct);
for i=1:length(fn)
    % load data files that exist. Rest will be empty
    fileName = [fn{i},'Name'];
    if any(strmatch(fileName,fn))
        try
            tmp = load(fullfile(dataStruct.dataFilePath,dataStruct.(fileName)));
            fnTmp = fieldnames(tmp);
            dataStruct.(fn{i}) = tmp.(fnTmp{1});
            
            % record named fields
            namedFieldsC = namedFieldsC + 1;
            namedFields{namedFieldsC,1} = fn{i};
            namedFields{namedFieldsC,2} = getVersion(dataStruct.(fileName));
            namesFields{namedFieldsC,3} = 1; % this is an indication of whether or not field is good to load
            namesFields{namedFieldsC,4} = dataStruct.(fileName);
            
            % get dependencies
            if isfield(dataStruct.(fn{i}),'dependencies')
                depListC = depListC + 1;
                dependencies{depListC,1} = fn{i};
                dependencies{depListC,2} = dataStruct.(fn{i}).dependencies;
                tempNames = fieldnames(dependencies{depListC,2});
                dependencies{depListC,3} = size(tempNames,1);
            end
        catch
            dataStruct.(fn{i}) = [];
        end
    end
end


%% version control

% if the number of named fields loaded is less than 3 then just end here
% because this is a new dataStruct
if namedFieldsC < 3
    return
end

% if no dependencies then this is an old analysis, so can end here
if depListC == 0
    warning('MAKILOADDATAFILE:NOVERSIONCONTROL',['--makiLoadDataFile: Warning, no version control in data file ' dataStruct.dataFileName]);
    return
end

% to check versions 1) start with list with the least dependencies and check
% version numbers of the 'Named' fields loaded in, 2) don't return 'Named' field with no version control

% list of numbers of deps
numDeps = NaN(size(dependencies,1),1);
for iDep = 1:size(dependencies,1)
    numDeps(iDep) = dependencies{iDep,3};
end
% min idx
[~,minIdx] = sort(numDeps);

% go over deps
depsFieldsC = 0;
for idx = minIdx'
    depsName = dependencies{idx,1};
    deps = dependencies{idx,2};
    
    depsNamePos = find(strcmp(depsName,namedFields(:,1)));
    
    depsFieldNames = fieldnames(deps);
    
    % go over deps fields
    
    for iFd = 1:length(depsFieldNames)
        
        % get recored version
        v = deps.(depsFieldNames{iFd});
        
        % add to list of dependency fields
        depsFieldsC = depsFieldsC + 1;
        depsFields{depsFieldsC} = depsFieldNames{iFd};
        
        % look up version in loaded table
        fieldPos = find(strcmp(depsFieldNames{iFd},namedFields(:,1)));
        
        % throw error if it's not there
        if isempty(fieldPos)
            error(['--makiLoadDataFile: Error: field ' depsName ' depends on field ' depsFieldNames{iFd} ', but it is not present']);
        end
        
        % see whether or not field is still good
        if ~namesFields{fieldPos,3}
            warning('MAKILOADDATAFILE:DEPENDENCYONBADFIELD',['--makiLoadDataFile: Warning, field ' depsName ' depends on field ' depsFieldNames{iFd} ', but this field has a version mismatch, will continue without loading ' depsName]);
            namesFields{depsNamePos,3} = 0;
            dataStruct.(depsName) = [];
        end
        
        % compare to current loaded version
        currentV = namedFields{fieldPos,2};
        
        % 3 possibilities, 1) version are equal, don't do anything, 2)
        % current version is higher, in this case lower the current version
        % number and attemps to load the correct file or 3) current version
        % is lower, then don't load the field
        if currentV > v
            warning('MAKILOADDATAFILE:CURRENTVERSIONTOOHIGH',['--makiLoadDataFile: Warning, field ' depsName ' depends on field ' depsFieldNames{iFd} ' version ' int2str(v) ', but the loaded version is ' int2str(currentV) '. Attemping to load version ' int2str(v)]);
            
            newFileName = changeVersion(namesFields{fieldPos,4},v);
            
            try
                tmp = load(fullfile(dataStruct.dataFilePath,newFileName));
                fnTmp = fieldnames(tmp);
                dataStruct.(fn{i}) = tmp.(fnTmp{1});
                namedFields{fieldPos,2} = v;
                namesFields{fieldPos,3} = 1; % this is an indication of whether or not field is good to load
                namesFields{fieldPos,4} = newFileName;
            catch
                warning('MAKILOADDATAFILE:COULDNOTLOADEARLIERVERSION',['--makiLoadDataFile: Warning, field ' depsName ' depends on field ' depsFieldNames{iFd} ' version ' int2str(v) ', but could not load that version, will continue without loading ' depsName]);
                namesFields{depsNamePos,3} = 0;
                dataStruct.(depsName) = [];
            end
            
            
        elseif currentV < v
            
            warning('MAKILOADDATAFILE:COULDNOTLOADEARLIERVERSION',['--makiLoadDataFile: Warning, field ' depsName ' depends on field ' depsFieldNames{iFd} ' version ' int2str(v) ', but could not load that version, will continue without loading ' depsName]);
            namesFields{depsNamePos,3} = 0;
            dataStruct.(depsName) = [];
            
        end
        
        
        
    end
    
    
end


% now check 'Named' filed againt those with version numbers, except for
% movieHeader, which won't have any dependencies AND the field with the
% most dependencies, idx is already set to the correct value

notToTest = {'movieHeader',depsName};

namesToTest = unique(namedFields(:,1));

namesToRemove = setdiff(namesToTest,depsFields);

for i = 1:length(notToTest)
    name = notToTest{i};
    namesNotToTest = strcmp(namesToRemove,name);
    namesToRemove = namesToRemove(~namesNotToTest);
end

% display warnings about fields that are to be removed
for i = 1:length(namesToRemove)
    warning('MAKILOADDATAFILE:FIELDHASNODEPENDIES',['--makiLoadDataFile: Warning, field ' namesToRemove{i} ' does not have any dependent fields, most likely bacuse of an analysis version mismatch, will not be loaded']);
    dataStruct.(namesToRemove{i}) = []; 
end
































