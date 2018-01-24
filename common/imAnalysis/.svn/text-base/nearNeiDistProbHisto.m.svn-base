function[pd,pdfilter]=nearNeiDistProbHisto(mpmstart);
% [pd,pdfilter]=nearNeiDistProbHisto(mpmstart);
% This function generates a probability distribution pd of nearest neighbor
% distances from the frame of an mpm file; as default, the distances are in
% steps of one pixel. In order to get a continuous probability (e.g. to be 
% used for MC simulations), the values are filtered to yield pdfilter
%
% INPUT: mpmstart = mpm file containing points
%
% OUTPUT: pd =  binned probability density (occurence probability) as a
%               function of distance in pixels
%         pdfilter =  the same, only filtered
%
% nearNeiDistProbHisto is ued by
%               - clusterQuantRipleyMC
%               - mpmMCaddsim
%
% Dinah Loerke, July 16, 2005
%


% remove zeros from mpm
mpm=[nonzeros(mpmstart(:,1)) nonzeros(mpmstart(:,2))];
[mpmx,mpmy]=size(mpm);

% calculate all distances
[dmat]=distanceMatrix(mpm,mpm);
nnvec = zeros(mpmx,1);

%extract nearest distances into vector nnvec
for n=1:mpmx
    pointvectemp=dmat(n,:);
    nnvec(n)=min(nonzeros(pointvectemp));
end

%some parameters for subsequent filtering are derived from the results
minval = min(nnvec);
meanval = mean(nnvec);
usemax = max(round(2*minval),round(1.4*meanval));

%histogram of values; center of bins chosen such that the final filtered
%vector and the unfiltered one both are 1:100 (which also implicitly
%limits the considered values of interpoint distances to 100 pixels!! 

shift=12;
%shift is chosen to somewhat match the value of sigma specified below
xvec = 1 : (100+shift-1);
ohist = hist(nnvec,xvec);
%exclude values which are too high, since they probably come from cells at
%the edge of the image where the true nearest neighbor is not seen
ohist(usemax:(100+shift-1)) = 0;
histsum = sum(ohist);
%probabilities are normalized occurence values
pdorig = ohist/histsum;
pd = pdorig(1:100);

% smooth this original result with Gaussian filter
xs=-shift:1:shift;
sig = 4;
amps=exp(-(xs.^2)/(2*(sig^2)));
namps=amps/sum(amps);
filtershape = namps;
%the above definition of the filter will introduce a shift to the filtered
%vector, this will is compensated by a counter-shift further below
[filtervec] = filter(filtershape,1,pdorig);
pdfilter = filtervec(shift:length(filtervec));


end % of function





function[m2]=distanceMatrix(c1,c2)
%this subfunction makes a neighbour-distance matrix for input matrix m1
%input: c1 (n1 x 2 points) and c2 (n2 x 2 points) matrices containing 
%the x,y coordinates of n1 or n2 points
%output: m2 (n1 x n2) matrix containing the distances of each point in c1 
%from each point in c2
[ncx1,ncy1]=size(c1);
[ncx2,ncy2]=size(c2);
m2=zeros(ncx1,ncx2);
for k=1:ncx1
    for n=1:ncx2
        d=sqrt((c1(k,1)-c2(n,1))^2+(c1(k,2)-c2(n,2))^2);
        m2(k,n)=d;
    end
end
end