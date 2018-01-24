function showMFrameTraj(MFT,color)
%showMFrameTraj Display the trajectories of speckles over multiframes.
%
% SYNOPSIS :
%    showMFrameTraj(MFT,color)
%
% INPUT :
%    MFT : An m-by-2*n matrix where 'm' is the number of starting points
%       of all the trajectories and 'n' is the number of frames. The '2*n'
%       columns have the form [y1 x1 ... yi xi ...] where (yi,xi) is the
%       position of the points in the i-th frame.
%    color : Matlab color : 'r', 'b', 'y' etc.

figure(gcf);

base = MFT(:,1:2);
plot(base(:,2),base(:,1),[color '.']);
for j = 1:2:size(MFT,2)-3
   dispV = (MFT(:,j+2:j+3)-MFT(:,j:j+1))*scale;
   quiver(base(:,2), base(:,1), dispV(:,2), dispV(:,1), 0, color);
   base = base+dispV;
end

hold off;

