
function visualizeClusters(receptorClusters,receptorPositions,plotScheme,circleSize,newFigure)

if newFigure
    figure
end
hold on

%get number of clusters
numClusters = size(receptorClusters,1);

%get cluster sizes and identify which cluster have more than one member
clusterSize = NaN(numClusters,1);
for iCluster = 1 : numClusters
    clusterSize(iCluster) = length(find(receptorClusters(iCluster,:)~=0));
end
clustersSizeNot1 = find(clusterSize~=1);

%plot clusters with one member
clusterMembers = receptorClusters(clusterSize==1,1);
switch plotScheme
    case 1
        plot(receptorPositions(clusterMembers,1),receptorPositions(clusterMembers,2),'LineStyle','none','Marker','o','Color',[0.5 0.5 0.5]);
    case 2
        plotcircle(receptorPositions(clusterMembers,:),circleSize*ones(length(clusterMembers),1),'Color',[0.5 0.5 0.5]);
end

%define color sequence for labeling clusters
colorSymbols = 'krgcmby';
maxClusterSize = max(clusterSize);
if maxClusterSize > 8
    colorSymbols(8:maxClusterSize) = 'y';
end

switch plotScheme

    case 1

        %go over all clusters of size greater than one and plot them
        for iCluster = clustersSizeNot1'
            clusterMembers = receptorClusters(iCluster,:);
            clusterMembers = clusterMembers(clusterMembers~=0);
            numMembers = clusterSize(iCluster);
            plot(receptorPositions(clusterMembers,1),receptorPositions(clusterMembers,2),'Marker','o','Color',colorSymbols(numMembers-1));
        end
        
    case 2
        
        for iClusterSize = 2 : maxClusterSize
            clusterMembers = receptorClusters(clusterSize==iClusterSize,1:iClusterSize);
            if ~isempty(clusterMembers)
                clusterMembers = clusterMembers(:);
                plotcircle(receptorPositions(clusterMembers,:),circleSize*ones(length(clusterMembers),1),'Color',colorSymbols(iClusterSize-1));
            end
        end
end
