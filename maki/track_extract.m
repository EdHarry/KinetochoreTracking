function [coord] = track_extract(sisterList,percent,maxDiffPercent)

%coord = struct('coord1',zeros(41,6),'coord2',zeros(41,6));
%coord = struct('coord1',0,'coord2',0);
coord=[];
clearCoord=0;
j=1;

for i=1:length(sisterList)
    
    coords1=sisterList(i).coords1;
    coords2=sisterList(i).coords2;
    
    if ~isempty(coords1)
        if (sum(~isnan(coords1(:,1))) >= (percent/100)*size(coords1,1) || sum(~isnan(coords2(:,1))) >= (percent/100)*size(coords2,1)) && abs((sum(~isnan(coords1(:,1)))./size(coords1,1)) - (sum(~isnan(coords2(:,1)))./size(coords2,1))).*100 <= maxDiffPercent
            
            if clearCoord==0
                clear coord
                clearCoord=1;
            end
           % i
            coord(j).coords1 = coords1;
            coord(j).coords2 = coords2;
            coord(j).idx = i;
            j=j+1;
        end
    end
end
end