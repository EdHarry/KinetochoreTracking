function [corrFacMatrix] = makeCorrFactorMatrix(imsiz, dist, samplesize, mask); 
% this function makes a matrix containing the correction factors for
% Ripley's circumference correction, for all points in an image of size
% imsiz, at specified sampling samplesize
% SYNOPSIS: [corrFacMatrix] = makeCorrFactorMatrix(imsiz, dist,samplesize, mask);
% 
% INPUT:    imsiz       = image size 
%           dist        = distance vector
%           samplesize  = size of the sampling, e.g. 10 pixel
%           mask (optional) = mask image of area of interest
%
% OUTPUT:   corrFacMatrix = has dimension (sx x sy x length(dist))
%           each (x,y,:) vector contains the correction factor for all
%           pertinent distances at the position x,y


%hw = waitbar(0,'calculate correction factor');

sx = imsiz(1);
sy = imsiz(2);
ld = length(dist);

sa = samplesize;
hsa = ceil(samplesize/2);

[xgrid,ygrid] = ndgrid(hsa:sa:sx,hsa:sa:sy); 

[lx,ly] = size(xgrid);

corrFacVector = nan*zeros(lx,ly,ld);
% initialize corrFacMatrix
corrFacMatrix = nan*zeros(sx,sy,ld);   

% if mask is entered, calculate euclidian distance map
if nargin==4
    edistmat = bwdist(1-mask);
end

for i=1:lx
    for k=1:ly
        px = xgrid(i,k);
        py = ygrid(i,k);
        
        % calculate correction factor (vector) for this x,y position
        if (nargin<4) 
            % if there is no mask, correction is based on rectangular
            % geometry; in this case, if the 'mirror' position already 
            % exists, use this value
            % min distance from x-edge
            edx = min(i,(lx+1-i)); 
            % min distance from x-edge
            edy = min(k,(ly+1-k));
            % since the indices i and k run from 1 to lx.ly, the first of
            % the four mirror values to be calculated will always be the
            % value at edx,edy    
            if (~isnan(corrFacVector(edx,edy,:)))
                corrFacVector(i,k,:) = corrFacVector(edx,edy,:);
            else
                corrFacVector(i,k,:) = circumferenceCorrFactor(px,py,dist,sx,sy);
            end
        % else if there's a mask, use it for calculation
        else
            % comment/uncomment between versions as desired
            
            % version 1 - is more precise, uses overlay of circle with mask
            % of area of interest, but is computationally slow
            % corrFacVector(i,k,:) = circumferenceCorrFactor_Mask(px,py,dist,sx,sy,mask);
            
            % version 2 - uses approximation of circle cut off by straight
            % line, at the distance corresponding the euclidian distance of
            % the center point form the edge of the area of interest; is
            % much faster
            % relevant euclidian distance
            ced = edistmat(px,py);
            corrFacVector(i,k,:) = circumferenceCorrFactor_MaskEU(px,py,dist,sx,sy,ced);
        end


        xlo = (i-1)*sa+1; xhi = min(sx,i*sa);
        ylo = (k-1)*sa+1; yhi = min(sy,k*sa);
        for r=1:length(dist)
            % fill matrix value with the corresponding vector position
            corrFacMatrix(xlo:xhi,ylo:yhi,r) = corrFacVector(i,k,r);
        end
        
     end % of for k
     
     %waitbar(i/lx);
     
end % of for i

%close(hw);


end % of function



%%======================================================================





function [m2]=DistanceMatrix(c1,c2);
%this subfunction makes a neighbour-distance matrix for input matrix c1
%(n1 x 2 points) and c2
%output: m2 (n1 x n1) matrix containing the distances of each point in c1 
%from each point in c2

[np1,sd1]=size(c1);
[np2,sd2]=size(c2);

m2=zeros(np1,np2);

for k = 1:np1
    for n = 1:np2
        d = sqrt((c1(k,1)-c2(n,1))^2+(c1(k,2)-c2(n,2))^2);
        m2(k,n)=d;
    end
end


end