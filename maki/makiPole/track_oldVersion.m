function track_oldVersion( jobVector, askDecon, basePath )
% EHarry Feb 2012

if nargin < 1 || strcmp(jobVector,'all')
    jobVector=[1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17];
end

if nargin < 2 || isempty(askDecon)
    askDecon = 1;
end

% if nargin < 3 || isempty(basePath)
%     basePath = uigetdir;
% end

if nargin < 3 || isempty(basePath)
    basePath = {};
    strResponse = '';
    while ~strcmp(strResponse,'done');
        basePath{end+1} = uigetdir;
        newbasePath = unique(basePath);
        if length(newbasePath) < length(basePath)
            warning(['Directory ',basePath{end},' has already been entered, ignoring.\n'])
        end
        basePath = newbasePath;
        strResponse = input('\nType "done" or press Enter to add another directory\n\n--> ', 's');
    end
end

if ~iscell(basePath)
    if ischar(basePath)
        basePath = {basePath};
    else
        error('Either enter a single directory as a string or multiple directories as a 1xn cell array of strings\n');
    end
end

fileList = [];

for i = 1:length(basePath)
    tempList = searchFiles('ome.tif',[],basePath{i},1);
    
    selectIdx = listSelectGUI(tempList(:,1),[],'move');
    % shorten tempList
    tempList = tempList(selectIdx,:);
    
    nJobs = length(tempList);
    
    if nJobs == 0
        %job = [];
        continue
    end
    
    
    
    if askDecon
        deconMovie(basePath{i},tempList);
    end
    
    fileList = [fileList;tempList];
end

% fileList = searchFiles('ome.tif',[],basePath,1);

%selectIdx = listSelectGUI(fileList(:,1),[],'move');
% shorten fileList
%fileList = fileList(selectIdx,:);

nJobs = size(fileList,1);

if nJobs == 0
    %job = [];
    return
end



job(1:nJobs) = struct('jobPath','/JobLogFiles',...
    'jobName','temp','dataStruct',[]);

%initialize progress display
progressText(0,'\nMoving Files into Place\n');

for iJob = 1:nJobs
    
    % for now: project name is rawMovieName minus extension
    rawMovieName = fileList{iJob,1};
    rawMoviePath = fileList{iJob,2};
    
    
    
    
    extIdx = regexp(rawMovieName,'ome.tif','split');
    extIdx = regexp(rawMovieName,extIdx{1},'end');
    
    projectName = rawMovieName(1:extIdx-1);
    
    % make folder for file and movie if necessary
    % For both test and hercules projects, rawMoviePath and
    % dataFilePath are identical. Thus, check whether the
    % rawMoviePath already contains the project name. If not, make
    % a subdirectory
    
    if any(findstr(rawMoviePath,projectName))
        dataFilePath = rawMoviePath;
        
    else
        dataFilePath = fullfile(rawMoviePath,projectName);
        mkdir(dataFilePath);
        rawMovieFile = searchFiles(projectName,'(log)|(_PRJ)|(hgsb)|(_history.txt)',rawMoviePath,0);
        if size(rawMovieFile,1) > 2
            errordlg(sprintf('%s is not a unique movie filename in the directory %s\nRename before continuing', ...
                rawMovieFile{1,1}, rawMovieFile{1,2}, ...
                'Non unique movie filename'));
            return;
        elseif size(rawMovieFile,1) == 2 % potentially have a raw movie and a deconvolved movie here
            numDecon = 0;
            for iList = 1:2
                if strcmp(rawMovieFile{iList,1}(end-7:end),'cmle.r3d')
                    numDecon = numDecon + 1;
                end
            end
            
            if numDecon == 1 % if exactly one decon movie was found then move both the files to the new dir
                for iList = 1:2
                    movefile(...
                        fullfile(rawMovieFile{iList,2},rawMovieFile{iList,1}),...
                        dataFilePath)
                end
                rawMoviePath = dataFilePath;
                %decon = 1;
            else
                errordlg(sprintf('%s is not a unique movie filename in the directory %s\nRename before continuing', ...
                    rawMovieFile{1,1}, rawMovieFile{1,2}, ...
                    'Non unique movie filename'));
            end
        else
            movefile(...
                fullfile(rawMovieFile{1,2},rawMovieFile{1,1}),...
                dataFilePath)
            rawMoviePath = dataFilePath;
            %decon = 0;
        end
        
    end
    
    job(iJob).dataStruct.rawMovieName = rawMovieName;
    job(iJob).dataStruct.rawMoviePath = rawMoviePath;
    job(iJob).dataStruct.dataFilePath = dataFilePath;
    
    progressText(iJob/(nJobs),sprintf('\nMoving %i files\n',nJobs));
    
end

job = makiMakeJobPlatformIndependent(job,'MCAINSH');

warning off
for iJob = 1:nJobs
    batchJobs(iJob).job = batch(['makiMakeJob_poleOME(7,2,job(',int2str(iJob),'),0,[],0);']);
end
warning on

progressText(0,'\nReading Movie Files\n');
for iJob = 1:nJobs
    batchJob = batchJobs(iJob).job;
    wait(batchJob);
    progressText(iJob/(nJobs),sprintf('\nRead %i Files\n',nJobs));
end

master_track_poleOME(job,jobVector,1,0);


end

