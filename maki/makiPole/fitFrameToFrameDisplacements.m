function bic = fitFrameToFrameDisplacements( data )
% EHarry October 2011
warning off
op = statset('MaxIter',1000);
obj1 = gmdistribution.fit(data,1,'Options',op);
obj2 = gmdistribution.fit(data,2,'Options',op);
warning on

mu = obj1.mu;
bic1 = obj1.BIC;
bic2 = obj2.BIC;

if bic1 > bic2 % if a 2 component model is better than a single then the data is not uniformly distributed
   bic = bic1 + 1e8; % in which case introduce a massive penilty 
else
    bic = bic1; % otherwise return the single model bic
end

if mu < 0.05 || mu > 0.3
    bic = bic + 1e10.*mu; % penilty for incorrect mu
end

end

