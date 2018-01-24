function [stimg,map]=RGStereo(i,j);
%RGSTEREO displays a red-green interlaced stereo picture from two input pictures
%
% SYNOPSIS rgStereo(i,j)
%
% INPUT    i,j  : gray value pictures of identical size
% OUTPUT   stimg: stereo picture
%          map  : colormap

%create a red/green colormap
cmap=[gray(128); gray(128)];
% 1...128 red
cmap(1:128,2:3)=0;
% 128..256 green
cmap(129:256,[1 3])=0;
%check for correct size
if(size(i) ~= size(j))
   return;
end;
% create a indexed image with values in correct range
i=double(gray2ind(i,128));
j=double(gray2ind(j,128));
% create the interlaced image
img(1:2:size(i,1),:)=i(1:2:size(i,1),:);
img(2:2:size(j,1),:)=j(2:2:size(j,1),:)+128;
iptsetpref('ImshowBorder', 'tight')
h=imshow(img,cmap);
stimg = img;
map=cmap;