function visuKinetochoreGraphs( subGraphs , coords )
%VISUKINETOCHOREGRAPHS visu kinetochore graphs
%   EHarry April 2012

ColOrd = get(gca,'ColorOrder');
[m,~] = size(ColOrd);

figure
plot3(coords(:,1),coords(:,2),coords(:,3),'b.');
hold on
for iSG = 1:length(subGraphs)
    ColRow = rem(iSG,m);
    if ColRow == 0
        ColRow = m;
    end
    Col = ColOrd(ColRow,:);
    spots = subGraphs(iSG).graph;
    coord = coords(spots,:);
    plot3(coord(:,1),coord(:,2),coord(:,3),'Color',Col);
end

end

