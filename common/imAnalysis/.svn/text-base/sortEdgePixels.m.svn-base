function edge_pixel_sorted = sortEdgePixels(edge_pixel)
%sortEdgePixels : Sort the pixels that define an edge.
%
% SYNOPSIS : edge_pixel_sorted = sortEdgePixels(edge_pixel)
%    Given the set of pixels that define an edge, return the same set of
%    pixels but sorted.
%
% INPUT :
%    edge_pixel : An m-by-2 matrix of the coordinates of the set of pixels
%       where m is the total number of pixels. The first column gives the x
%       (horizontal) coordinate and the second column gives the y (vertical)
%       coordinate.

if ~isnumeric(edge_pixel) | size(edge_pixel,2) ~= 2
   error('The input should be an m-by-2 matrix.');
end

%We choose the leftmost pixel as the first pixel to start. If it is not
% unique, the choose among them the upper-most.
[hmin,ind] = min(edge_pixel(:,1));
ind = find(edge_pixel(:,1) == hmin);
[vmin,ind] = min(edge_pixel(ind,2));

edge_pixel_sorted = edge_pixel(ind,:);
edge_pixel(ind,:) = [];

i=1;
i_s=1;
found =1;
while size(edge_pixel,1) > 1 && found ==1
   %find first 4 connected neighbour
   while i <= size(edge_pixel,1) && ...
      sqrt((round(edge_pixel(i,1))-round(edge_pixel_sorted(i_s,1)))^2 + ...
      (round(edge_pixel(i,2))-round(edge_pixel_sorted(i_s,2)))^2) > 1
      i=i+1;
   end
   if i>size(edge_pixel,1);
      %now try to find a 8 connected neighbour
      i=1;
      while  i <= size(edge_pixel,1) && ...
         sqrt((round(edge_pixel(i,1))-round(edge_pixel_sorted(i_s,1)))^2 + ...
         (round(edge_pixel(i,2))-round(edge_pixel_sorted(i_s,2)))^2) > sqrt(2)
         i=i+1;
      end
      if i>size(edge_pixel,1);
         %no next pixel was found
         found=0;
      end
   end
   if found~=0;
      i_s=i_s+1;
      edge_pixel_sorted(i_s,:) = edge_pixel(i,:);
      edge_pixel(i,:)=[];
      i=1;
   end
end
