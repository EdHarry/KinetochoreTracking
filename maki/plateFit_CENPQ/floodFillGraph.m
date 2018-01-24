function  subGraphs  = floodFillGraph( c )
% EHarry Sept 2011

numS = size(c,1);
idx = 1:numS;

col = 1;

loop=1;
while loop
    colTemp=[];
    for q = 1:length(col)
        connected = find(c(:,col(q)));
        colTemp(end+1:end+length(connected)) = connected; % inital run
    end
    col2 = unique([col,colTemp]);
    if length(col2) == length(col)
        loop=0;
    end
    col = col2;
end

subGraphs(1).graph = col;

notIncluded = idx(~ismember(idx,col));

while ~isempty(notIncluded)
    col = notIncluded(1);
    
    loop=1;
    while loop
        colTemp=[];
        for q = 1:length(col)
            connected = find(c(:,col(q)));
            colTemp(end+1:end+length(connected)) = connected; % inital run
        end
        col2 = unique([col,colTemp]);
        if length(col2) == length(col)
            loop=0;
        end
        col = col2;
    end
    
    subGraphs(end+1).graph = col;
    
    colAll=[];
    for w = 1:length(subGraphs)
        colTemp = subGraphs(w).graph;
        colAll(end+1:end+length(colTemp)) = colTemp;
    end
    
    notIncluded = idx(~ismember(idx,colAll));
end

end

