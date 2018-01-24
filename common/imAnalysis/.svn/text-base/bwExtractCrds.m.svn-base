function [xi,yi] = bwExtractCrds(bw)
%BWEXTRACTCRDS extracts the non zero pixel coordinates of a binary mask
%
% SYNOPSIS [xi,yi] = bwExtractCrds(bw)
%
% INPUT bw : binary image with values {0,1}
%
% OUTPUT xi : vector with x_1 coordinates
%        yi : vector with x_2 coordinates

nPts = sum(sum(bw));
xi = zeros(1,nPts);
yi = xi;

n = 0;

for(i=1:size(bw,1))
   for(j=1:size(bw,2))
      if(bw(i,j))
         n = n+1;
         xi(n) = j;
         yi(n) = i;
      end;
   end;
end;
