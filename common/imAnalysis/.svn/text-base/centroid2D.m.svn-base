function ce = centroid2D(img)
% CENTROID compute the centoid of a gray value patch
%
% SYNOPSIS ce = centroid2D(img)
%
% INPUT img : an image patch matrix
% 
% OUTPUT ce : vector with the centroid coordinates

% c 19/04/00

s = size(img);
cx = 0;
cy = 0;
for l = 1 : s(2);
   cx = cx + sum(img(:,l).*l);
end
for l = 1 : s(1);
   cy = cy + sum(img(l,:).*l);
end
sTot=sum(img(:));
ce=[cx cy]/sTot;