function [corfac] = circumferenceCorrFactor_Mask(xx,yy,rr,msx,msy,mask)
% circumference correction calculates a vector containing the correction factor
% (for edge correction in Ripley's k-function) for values of rr
% circumference correction: fraction of circumference of circle centered at
% point P=(xx,yy) with radius rr (inside rectangular image) falling into the 
% area of interest defined by the mask image 
%
% Notes: - this function assumes that rr is a vector
%        - this version of the function allows rr to be larger than ms/2

% SYNOPSIS
% [corfac]=circumferenceCorrFactor_Mask(xx,yy,rr,msx,msy,mask)
% INPUT:    xx      x-position 
%           yy      y-position
%           rr      distance (vector)
%           msx     image size in x
%           msy     image size in y
%           mask    mask of area of interest, if it is not the rectangular
%                   image
%
% OUPUT:    corfac  correction factor; vector of same length as rr 
%       
% Dinah Loerke, January 29, 2008

% the dimensions of the mask have to be the same as msx, msy
[mx,my] = size(mask);
if (mx~=msx) | (my~=msy)
    error('image and mask size don''t agree');
end

% due to the grid-based/sample approach to the correction factor problem,
% it can happen occasionally that the chosen coordinates xx,yy lie OUTSIDE
% the area mask; there are two possibilities for this case, on the one
% hand, cofac values could be set to nan, on the other hand, cofac values
% could be set to the values for the CLOSEST point to xx,yy inside the area

px = round(xx);
py = round(yy);

if mask(px,py)==0
    inpoints = find(mask);
    [inpx,inpy] = ind2sub([msx,msy],inpoints);
    dist = DistanceMatrix([xx,yy],[inpx inpy]);
    pdist = find(dist==min(dist));
    xx = inpx(pdist);
    yy = inpy(pdist);
end
    

rl = length(rr);
%correction factor is initialized
corfac=ones(rl,1);
% maximum radius = rs
rs = max(rr);

% make stack of circle masks for appropriate rr-distances
[circleImX, circleImY] = ndgrid(-rs:rs,-rs:rs);
circleDist = sqrt(circleImX.^2 + circleImY.^2);
% stack of masks which contain the ring between distances rr(s) and rr(s+1)
rr2 = [0,rr];
for s=1:rl
    maskStack(:,:,s) = ((circleDist<rr(s)) & (circleDist>=rr2(s)));
end

% for the mask image, the x/y coordinates are NOT switched here
% pad the image at the edges with rs
maskWE = zeros(msx+2*rs,msy+2*rs);
maskWE(rs+1:rs+msx,rs+1:rs+msy) = mask;

% loop over all distance values
for ri = 1:rl
    
    % piece of the analysis area cut around with rs around the current
    % point of interest
    % for the mask image, the x/y corrdinates are switched
    currentMaskCutout = maskWE(xx:xx+2*rs,yy:yy+2*rs);

    % product
    circleMask = maskStack(:,:,ri);
    prodIm = currentMaskCutout.*circleMask;
    % overlap = number of maskarea pixels found within this circle
    % radius
    mapix_num = max(sum(prodIm(:)),1);
    % total circle area for this distance
    %ta = pi*ri^2;
    ta = sum(circleMask(:));
    % correction factor is fraction of analysis area inside the circle
    corfac(ri)=mapix_num/ta;

end % of for ri
    

end  % of function


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
