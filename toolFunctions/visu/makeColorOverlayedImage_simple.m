% function [ spotImage, colorImage, map ] = makeColorOverlayedImage_simple( coords, psfSigma, sizeX , sizeY, sizeZ )
function [spotImage, colorImage, map] = makeColorOverlayedImage_simple( coords, psfSigma, sizeX , sizeY, sizeZ )
%MAKECOLOROVERLAYEDIMAGE Summary of this function goes here
%   Detailed explanation goes here

[~,spotImage] = placeGausianSpots( coords , psfSigma, [sizeX sizeY sizeZ]);

map = lines(size(coords,1));

colorImage = zeros(sizeX,sizeY,sizeZ,3);
colorImage(:,:,:,2) = 1; 

%colorImage = zeros(sizeX,sizeY,sizeZ,3);
% for i = 1:size(coords,1)
%     imageCoords = placeGausianSpots( coords(i,:) , psfSigma, [sizeX sizeY sizeZ]);
%     colorImage(imageCoords(:,1),imageCoords(:,2),imageCoords(:,3),1) = colorImage(imageCoords(:,1),imageCoords(:,2),imageCoords(:,3),1) + map(i,1);
%     colorImage(imageCoords(:,1),imageCoords(:,2),imageCoords(:,3),2) = colorImage(imageCoords(:,1),imageCoords(:,2),imageCoords(:,3),2) + map(i,2);
%     colorImage(imageCoords(:,1),imageCoords(:,2),imageCoords(:,3),3) = colorImage(imageCoords(:,1),imageCoords(:,2),imageCoords(:,3),3) + map(i,3);
% end
% for i = 1:size(coords,1)
%     imageCoords = placeGausianSpots( coords(i,:) , psfSigma, [sizeX sizeY sizeZ]);
%     ind = sub2ind([sizeX sizeY sizeZ],imageCoords(:,1),imageCoords(:,2),imageCoords(:,3));
%     temp = colorImage(:,:,:,1);
%     temp(ind) = temp(ind) + map(i,1);
%     colorImage(:,:,:,1) = temp;
%     temp = colorImage(:,:,:,2);
%     temp(ind) = temp(ind) + map(i,2);
%     colorImage(:,:,:,2) = temp;
%     temp = colorImage(:,:,:,3);
%     temp(ind) = temp(ind) + map(i,3);
%     colorImage(:,:,:,3) = temp;
% end

end

