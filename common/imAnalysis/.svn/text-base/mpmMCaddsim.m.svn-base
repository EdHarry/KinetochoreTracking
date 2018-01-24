function[mpm]=mpmMCaddsim(mpmstart,imsize, numcells,pd)
%[mpm]=mpmMCaddsim(mpmstart,imsize, numcells,pd)
%This function uses MonteCarlo simulation to add cells to an existing mpm 
%of a cell distibution. The existing mpm is used to calculate the probability 
%distribution of the distance of the nearest neighbor; new cells are placed
%into the image accordingly
%
%INPUT: mpmstart = the original distribution of cells
%       imsize = [xsize ysize] vector of image size
%       numcells = number of cells to be added to the image
%       pd = probability density; if this is not entered, it's calculated
%       from mpmstart

%OUTPUT: mpm = the new mpm with added cells
%
% strategy of the function : unless provided, determine nearest neighbor 
% distance probability from the input mpm (histogram of actual nearest 
% neighbor distances); in the next step, place specified
% number of new points at locations defined by these probabilities
%
% mpmMCaddsim is used by clusterQuantRipleyMC and variations thereof
%             uses nearNeiDistProbHisto
%
% Dinah Loerke, July 16, 2005


%=========================================================================
% step 1: calculate nearest neighbor probability from frame of mpm unless 
% pd is already specified as an input
%=========================================================================
if nargin<4
    [pd,pdfilter]=nearNeiDistProbHisto(mpmstart);
else
    pdfilter = pd;
end

%pd is probability as a function of d vector (in steps of 1 pixel)

mpm=[nonzeros(mpmstart(:,1)) nonzeros(mpmstart(:,2))];

%loop over number of points that are supposed to be added to image
for i=1:numcells
    
    %celldist contains distance from nearest cell rad point
   
    c=0;
    %while the generated point does not yet fulfill the requirements
    while (c==0)
        %generate random point
        pointx = ceil(imsize(1)*rand(1)) ;
        pointy = ceil(imsize(2)*rand(1)) ;
        %calculate nearest neighbor distance to existing cells
        [currCellDist]=distanceMatrix([pointx pointy],mpm);
        currNNdist=min(nonzeros(currCellDist));
        
        %determine probability for this distance (if it's in the right
        %ballpark, i.e. no more than 100 pix)
        if (currNNdist < 100)
            %occurence probability for this particular distance
            currProb = pdfilter(ceil(currNNdist));
            %generate random variable between 0 and 1
            randomVar = rand(1);
            %do MC test
            if ( randomVar <  currProb )
                c=1;
            end
        end % of if
        
    end % of while
    
    newpointvec = [pointx pointy];
    mpm=[mpm; newpointvec];

    
end % of for

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


