function data=pushstamp3d(data,patch,center,symmetric,add)
%PUSHSTAMP3D writes a 3Dsubimage into a larger image 
%
% SYNOPSIS data=pushstamp3d(data,patchSize,center)
%
% INPUT data   : 3D data
%            patchSize  : size of patch
%            center     : pixel coords of center of patch in 3D data
%            symmetric  : (opt) if part of the patch would fall outside of the 
%                           img: whether to cut accordingly on the other side 
%                           [{0}/1] 
%            add        : (opt) whether to replace or to add patch to image
%                           or whether to subtract [{0}/1/-1]
%            
%
% OUTPUT data : 3D data 

% c: 18/6/01	dT

% test input, assign defaults
if nargin < 4 || isempty(symmetric)
    symmetric = 0;
end

if nargin < 5 || isempty(add)
    add = 0;
end

% find patch in data
ds=size(data);
patchSize = size(patch);
if length(patchSize) == 2
    patchSize(3) = 1;
end
hl = patchSize/2 - 0.5;

% find extension towards 0
hx1=min([center(1)-1,hl(1)]);
hy1=min([center(2)-1,hl(2)]);
hz1=min([center(3)-1,hl(3)]);

% find extension towards inf
hx2 = min([hl(1),ds(1)-center(1)]);
hy2 = min([hl(2),ds(2)-center(2)]);
hz2 = min([hl(3),ds(3)-center(3)]);

% make symmetric, if necessary
if symmetric
   [hx1,hx2] = deal(min(hx1,hx2));
   [hy1,hy2] = deal(min(hy1,hy2));
   [hz1,hz2] = deal(min(hz1,hz2));
end

switch add
    case 1 % add
        data(center(1)-hx1:center(1)+hx2,...
            center(2)-hy1:center(2)+hy2,...
            center(3)-hz1:center(3)+hz2) =...
            data(center(1)-hx1:center(1)+hx2,...
            center(2)-hy1:center(2)+hy2,...
            center(3)-hz1:center(3)+hz2) + ...
            patch(hl(1)+1-hx1:hl(1)+1+hx2,...
            hl(2)+1-hy1:hl(2)+1+hy2,...
            hl(3)+1-hz1:hl(3)+1+hz2);
    case -1 % subtract
        data(center(1)-hx1:center(1)+hx2,...
            center(2)-hy1:center(2)+hy2,...
            center(3)-hz1:center(3)+hz2) =...
            data(center(1)-hx1:center(1)+hx2,...
            center(2)-hy1:center(2)+hy2,...
            center(3)-hz1:center(3)+hz2) - ...
            patch(hl(1)+1-hx1:hl(1)+1+hx2,...
            hl(2)+1-hy1:hl(2)+1+hy2,...
            hl(3)+1-hz1:hl(3)+1+hz2);
    otherwise % replace
        data(center(1)-hx1:center(1)+hx2,...
            center(2)-hy1:center(2)+hy2,...
            center(3)-hz1:center(3)+hz2) =...
            patch(hl(1)+1-hx1:hl(1)+1+hx2,...
            hl(2)+1-hy1:hl(2)+1+hy2,...
            hl(3)+1-hz1:hl(3)+1+hz2);
end
