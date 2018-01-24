function  [mov, map] = stack2Movie(stack,axesH,nCol,nRep,fps)
%STACK2MOVIE converts a stack of intensity images to a MATLAB movie 
% and displays it nRep times on the axis given by axesH with fps
% frames per second
%
% SYNOPSIS [mov, map] = stack2Movie(stack,axesH,nCol,nRep,fps)
%
% INPUT stack : image stack
%       axesH: axes handle
%       nCol : number of entries in the colormap
%       nRep : (optional) number of repetitions (default: 1)
%       fps : (optional) frames per second (default: 12)
%
% OUTPUT mov : movie data
%        map : colormap used
%
% SEE ALSO movie

if(nargin < 5)
   fps = 12;
end;

if(nargin < 4)
   nRep = 1;
end;

axes(axesH);

[stackHeight, stackWidth, stackDepth] = size(stack);

[indImg, map] = gray2ind(stack(:,:,1),nCol);

for iFm = 2:stackDepth
   indImg = cat(4,indImg,gray2ind(stack(:,:,iFm),nCol));
end;
mov = immovie(indImg,map);

if nRep > 1
   movie(mov,nRep-1,fps)
end;


