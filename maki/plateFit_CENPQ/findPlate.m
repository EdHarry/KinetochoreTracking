function subGraphs = findPlate( coord,cutOff,minNum )
% EHarry Sept 2011

subGraphsTemp  = floodFillGraph( makeGraph( coord,cutOff ) );
subGraphs = struct('graph',[]);

temp=[];
idx=1;
for i = 1:length(subGraphsTemp)
    if length(subGraphsTemp(i).graph) > minNum-1
        subGraphs(idx).graph = subGraphsTemp(i).graph;
        temp(idx) = length(subGraphs(idx).graph);
        idx = idx + 1;
    end
end

if ~isempty(temp)
    [~,i] = sort(temp,2,'descend');
    subGraphs = subGraphs(i); % sort the reults in decending order
end

end

