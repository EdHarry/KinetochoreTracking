function c = makeGraph( coords,cutOff )
%MAKEGRAPH creates connectivity graph of spots defined in coords
%   INPUT
%   coords: (n by 3) matrix of spot coords
%   cutOff: cutOff distance over which a connection is not allowed
%
%   OUTPUT
%   c: connectivity matrix of spots
% EHarry Sept 2011

numSpots = size(coords,1);

c = zeros(numSpots); % set inital matrix to completely disconected

d = createSparseDistanceMatrix(coords,coords,cutOff); % distance matrix, zeros after cutoff
d = full(d);
d(d==0) = 1000; % set zeros to 1000, a big number

for j = 1:numSpots % set all self-distance to 1000
    d(j,j) = 1000;
end


numCon = 1;
while numCon > 0
    numCon = 0;
    [~,i] = min(d,[],2); % min indicies
    
    for j = 1:numSpots
        if d(j,i(j)) < 1000
            numCon = numCon + 1;
            c(j,i(j)) = 1; % set indicies to 1 in contact matrix, but only if the distance is less than 1000 (not one that has been excluded)
        end
    end
    
    for j = 1:numSpots
        d(j,i(j)) = 1000; % set found indicies to 1000
    end
end

end

