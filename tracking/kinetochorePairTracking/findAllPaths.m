function paths = findAllPaths( conflicts )
%FINDALLPATHS Summary of this function goes here
%   Detailed explanation goes here


n = length(conflicts);
m = 0;
for i = 1:n
    mm = max(conflicts{i});
    if mm > m
        m = mm;
    end
end

paths = {};
for i = 1:m
    paths = [paths; {i}];%#ok<AGROW>
end


for time = 2:n+1
    l = length(paths);
    for i = 1:l
        newPath = [paths{i} paths{i}(end)];
        paths = [paths; newPath];%#ok<AGROW>
    end
    p = conflicts{time-1};
    for i = 1:length(p)-1
        for j = i+1:length(p)
            states = [];
            for k = 1:l
                states = [states; paths{k}(end)];%#ok<AGROW>
            end
            t1 = find(states==p(i));
            t2 = find(states==p(j));
            for k = 1:length(t1)
                newPath = [paths{t1(k)} p(j)];
                paths = [paths; newPath];%#ok<AGROW>
            end
            for k = 1:length(t2)
                newPath = [paths{t2(k)} p(i)];
                paths = [paths; newPath];%#ok<AGROW>
            end
        end
    end
    
    paths(1:l) = [];
end

paths(1:m) = [];

end

