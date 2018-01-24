function f = meanNNDist(M);
%MEANNNDIST computes the mean nearest neighbour distance of an N-dim point set
%
% SYNOPSIS mean = meanNNDist(M)
%
% INPUT:  M: an S x T matrix where the # of columns is the # of dimensions
%            and  the # of rows is the # of points
%
% c: 20/9/99	dT

% compute the distance matrix
[rows,cols]=size(M);
for i = 1:rows
   for j = 1:rows
      if(j<i)
         DIST(i,j) = norm(M(j,:)-M(i,:));
      end;
      if(j>i)
         DIST(i,j-1) = norm(M(j,:)-M(i,:));
      end;
   end;
end;
% minV is a vector of the minimum element of each row of DIST
minV=min(DIST,[],2);
% compute the mean value over minV
f = mean(minV);
 