function [poleCoords,sList] = poleCoords
%Wrapper for converting to pole coord system
% EHarry March 2011

poleCoords=0; % placeholder values for poleCoords and sList if script returns no results
sList=0;
clearPoleCoords=0; % flag for clearing poleCoords and sList

basePath = uigetdir('','Please select directory of movie');   % basepath to movie

fileListPlaneFit = searchFiles('planeFit.*mat','',basePath,1,'new');  % search for planeFit mat files

counter=1; % struct counter
for i=1:length(fileListPlaneFit)
    
    fileListSisterList = searchFiles('sisterList.*mat','',fileListPlaneFit{i,2},1,'new');
    fileListTracks = searchFiles('tracks.*mat','',fileListPlaneFit{i,2},1,'new');
    fileListPole = searchFiles('poleCoord.*mat','',fileListPlaneFit{i,2},1,'new'); % try to load sisterLists, tracks and poleCoords for each subdirectory with a planeFit
    
    eflag=0;
    
    if ~isempty(fileListPole) && ~isempty(fileListSisterList)
        load(fullfile(fileListPole{2},fileListPole{1})); % load the poleCoords from file if it exists
        load(fullfile(fileListSisterList{2},fileListSisterList{1})); % load sisterList as well
    else % else try and load the sisterList and tracks for each planeFit and run the conversion
        
        if ~isempty(fileListSisterList) && ~isempty(fileListTracks) % only continue if there is a sisterList and a tracks for each found planeFit
            load(fullfile(fileListPlaneFit{i,2},fileListPlaneFit{i,1})); % load planeFit
            load(fullfile(fileListSisterList{2},fileListSisterList{1})); % load sisterList
            load(fullfile(fileListTracks{2},fileListTracks{1}));  % load tracks
            [eflag,poleCoord] = cylindericalTransform(planeFit,tracks,sisterList); % run conversion
        else
            continue % if no sisterList or tracks for this planeFit then skip to the next loop iteration
        end
    end
    
    
    if ~eflag % only try and save/return results if no errors in conversion (or if poleCoords were loaded from file)
        
        if ~clearPoleCoords % clear poleCoords and sList to save results as structs
            clear poleCoords
            clear sList
            clearPoleCoords=1; % set flag so only one clear
        end
        
        if isempty(fileListPole) % save poleCoord results to poleCoord.mat file if the file doesn't exist
            save(fullfile(fileListPlaneFit{i,2},['poleCoord',fileListPlaneFit{i,1}(9:end)]),'poleCoord'); % save results to a .mat file in the directory
        end
        
        poleCoords(counter) = poleCoord; % save results to output
        sList(counter).sisterList = sisterList;
        counter=counter+1; % inc counter
    end
end
end


