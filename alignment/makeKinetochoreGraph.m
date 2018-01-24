function subGraphs = makeKinetochoreGraph( coords , verbose)
%MAKEKINETOCHOREGRAPH makes connectivity graph (or connectivity matrix) between features in a movie
%and detects disconected subgraphs
% EHarry April 2012

warning('off','ROBUSTMEAN:INSUFFICIENTDATA') % turn the warning off since it'll be used a lot, if there is insufficient data the oulier list will be empty anyway so it's fine

if nargin < 2
    verbose = 0;
end

numCoords = size(coords,1);

% if less than 4 coords then do a simple return
if numCoords < 4
    subGraphs(1:numCoords) = struct('graph',[]);
    for i = 1:numCoords
        subGraphs(i).graph = i;
    end
    return
end

% make an internal distance matrix
disMat = createDistanceMatrix(coords,coords);

% sort disMat
[~,sortIdx] = sort(disMat,2);

% get nearst neighbours
% nnMat = sortMat(:,2);
nnIdx = sortIdx(:,2);

% pair up nn
nnList = [(1:numCoords)',nnIdx];
nnList = sort(nnList,2);

% get unique pairs
nnList = unique(nnList,'rows');

% get the list of these nn distances
nnDis = diag(disMat(nnList(:,1),nnList(:,2)));

% detect outliers
[~,~,~,outliers] = robustMean(nnDis);

% remove outliers
nnList(outliers,:) = [];
nnDis(outliers) = [];

% remove any more than a sutable maxDis, say 3 microns
maxDis = 4;
nnList(nnDis>maxDis,:) = [];

% get the number of nn
numNN = size(nnList,1);

% make connectivity matrix between these nn
tempGraph = zeros(numCoords);
for i = 1:numNN
    tempGraph(nnList(i,1),nnList(i,2)) = 1;
    tempGraph(nnList(i,2),nnList(i,1)) = 1;
end

% get the connected subGraphs
subGraphs  = floodFillGraph( tempGraph );

% now loop and join up subGraphs that have acceptably small inter-graph
% connections
% nn is the idx to the nn of spots to join sub-graphs
nn = 2;
% while nn < floor(numCoords/6)
while nn < 3
    if 0
        close all
        figure,hold on
        for iSG = 1:length(subGraphs)
            spots = subGraphs(iSG).graph;
            coord = coords(spots,:);
            plot3(coord(:,1),coord(:,2),coord(:,3),'r-');
        end
        pause(0.5)
    end
    % increase nn
    nn = nn + 1;
    
    % get the size of all the subgraphs, don't process graphs smaller than
    % 5 spots (let larger graphs connect to them) UNLESS there are none bigger than that
    % since the subgraphs are sorted by size only need to test the first
    % one
    minGraphSize = 5;
    sizeLargest = length(subGraphs(1).graph);
    
    if sizeLargest < minGraphSize
        sizeRestriction = 0;
    else
        sizeRestriction = 1;
    end
    
    % go over subGraphs
    for iSubGraph = 1:size(subGraphs,1)
        % flag for whether or not new connections were made
        addedNew = false;
        
        spots = subGraphs(iSubGraph).graph;
        
        % test for size restriction, if these is and this graph is too
        % small skip it
        if sizeRestriction && length(spots) < minGraphSize
            continue
        end
        
        % get the nn of the spots
        nnIdx = (sortIdx(spots',nn-1))';
        
        % get the current distances for these spots
        currentDis = (diag(disMat(nnIdx,spots)))';
        
        %         % get the intra-subgraph internal distances
        %         intraDis = disMat(spots,spots);
        %         intraDis = triu(intraDis,1);
        %
        %         currentDis = (intraDis(intraDis ~= 0))';
        
        % get next nn idxs to these spots
        nnIdx = (sortIdx(spots',nn))';
        
        % get any of these nn that aren't already in the subGraph
        nnIdxNew = setdiff(nnIdx,spots);
        
        % if there are any then potentially add these connections to the
        % tempGraph if the new connections fall withing the scatter of the
        % current inter spot distances in the subgraph
        for idxNew = nnIdxNew
            % get the new idx pair
            newPair = [idxNew spots(nnIdx==idxNew)];
            
            % get the distance for this connection
            dis = disMat(newPair(1),newPair(2));
            
            % use robust-mean
            [meanDis,stdDis] = robustMean(currentDis);
            
            if verbose
                disp([int2str(iSubGraph) ' ' num2str(meanDis) ' ' num2str(stdDis) ' ' num2str(dis)])
            end
            
            % add connectivity to matrix if the new dis within a std of the
            % mean
            if abs(dis-meanDis) < 3*stdDis && dis < maxDis
                tempGraph(newPair(1),newPair(2)) = 1;
                tempGraph(newPair(2),newPair(1)) = 1;
                % flag that a new connection was make
                addedNew = true;
                % update the current list of distances
                currentDis = [currentDis dis];
            end
        end
        
        % if a new connection was make then make a new subGraph list, break
        % and restart the loop
        
        if addedNew
            subGraphs  = floodFillGraph( tempGraph );
            nn = nn - 1; % decrease nn so it remains the same for the next loop
            break
        end
        
    end
    if verbose
        disp(['num subgraphs: ',int2str(length(subGraphs)) ', nn: ' int2str(nn)])
    end
    
end


% now get sizes of the subgraphs
sizes = NaN(length(subGraphs),1);
for iSG = 1:length(subGraphs)
    sizes(iSG) = length(subGraphs(iSG).graph);
end

% if any size is greater than one, continue
if ~any(sizes > 1)
    return
end

% cluser spots in subgraphs by mean coord
newCoords = NaN(length(subGraphs),3);
newCoord2OrigMap = cell(length(subGraphs),1);
for iSG = 1:length(subGraphs)
    spots = subGraphs(iSG).graph;
    coord = coords(spots,1:3);
    % use robust mean
    newCoords(iSG,1) = robustMean(coord(:,1));
    newCoords(iSG,2) = robustMean(coord(:,2));
    newCoords(iSG,3) = robustMean(coord(:,3));
    % fill in map
    newCoord2OrigMap{iSG} = spots;
end

% rerun makeKinetochoreGraph on new coords
subGraphs = makeKinetochoreGraph( newCoords,verbose );

% substitute in original coords
for iSG = 1:length(subGraphs)
    spots = subGraphs(iSG).graph;
    originalSpots = [];
    for spot = spots
        mappedOriginalSpots = newCoord2OrigMap{spot};
        originalSpots(end+1:end+length(mappedOriginalSpots)) = mappedOriginalSpots;
    end
    subGraphs(iSG).graph = originalSpots;
end

warning('on','ROBUSTMEAN:INSUFFICIENTDATA')
end

