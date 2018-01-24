function trajDataOverlyResults( trajData, showFullTracks )
%TRAJDATAOVERLYRESULTS overlay trajData results onto of image
%   EHarry Oct 2012

if nargin < 2
    showFullTracks = 0;
end

for iCell = 1:length(trajData)
    % select single cell
    traj = trajData(iCell);
    
    % display cell info
    disp(traj.cellName)
    
    % load corresponding dataStruct
    dataStruct = makiLoadDataFile('MCAINSH',fullfile(traj.cellPath,traj.cellName));
    
    % get indexes of extracted kinetochores
    sisIdx = catStruct(1,'traj.sisterList.idx');
    
    % delete appropriate spots if not showing full tracks
    if ~showFullTracks
        for iSis = 1:length(sisIdx)
            toDel = isnan(traj.sisterList(iSis).coords1(:,1)) | isnan(traj.sisterList(iSis).coords2(:,1));
            dataStruct.sisterList(sisIdx(iSis)).coords1(toDel,:) = NaN;
            dataStruct.sisterList(sisIdx(iSis)).coords2(toDel,:) = NaN;
        end
    end
    
    % show the image
    overlaySistersOnImage( dataStruct, sisIdx );
    
    % wait for user input
    strResponse = '';
    while ~strcmp(strResponse,'done')
        strResponse = input('\nEnter "done" to load the next cell\n\n--> ', 's');
    end
end

end

