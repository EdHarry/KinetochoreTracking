function [img, dummy] = ppMonoFgGrab(cam,dummy)
%PPMONOFGGRAB grabs an image with the 2xPicPortMono framegrabber
%
% SYNOPSIS  [img, fRoi] = ppMonoFgGrab(roi,rgb)
%
% INPUT     cam: integer (0=left camera, 1=right camera)
% OUTPUT    img: intensity image (m x n matrix)
%
% NOTE: dummy arguments just for compatibility

dummy=[];
[status, result] = mexFgInterface('grab');

%use left cam only
if((cam==0) | (cam == 1))
   img= squeeze(result(:,:,cam+1))';
end;
return;  