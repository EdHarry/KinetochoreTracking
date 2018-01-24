function nbIdx =neighbours(pSet1,pSet2,nbDist);
% NEIGHBOURS identifies neighbours of two point sets
%
% SYNOPSIS : [nbIdx]=neighbours(pSet1,pSet2,nbDist)
%
% INPUT  pSet1,pSet2: two (n X m matrix) point sets with n=dim of a point
%                     and m = number of points
%        nbDist     : scalar neighbour distance
%
% OUTPUT nbIdx      : (2 x s) matrix, first row indexes into pSet1
%                     second row indexes into pSet2, each column pair
%                     corresponds to two neighbours

lSet1=size(pSet1,2);
lSet2=size(pSet2,2);
nbIdx=[0 ;0];
for i=1:lSet1
   % for each point of pSet1 find neighbours in pSet2
   distVec=zeros(1,lSet2);
   for d=1:size(pSet1,1)
      distVec=distVec+(pSet2(d,:)-pSet1(d,i)).^2;
   end;
   lenDist=sqrt(distVec)<nbDist;
   idPs2=find(lenDist);
   nbIdx=[nbIdx(1,:) i*ones(1,length(idPs2)); nbIdx(2,:) idPs2];
end;
if(size(nbIdx,2)>1)
   nbIdx=[nbIdx(1,2:length(nbIdx));nbIdx(2,2:length(nbIdx))];
else
   nbIdx=[];
end;


