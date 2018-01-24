function [bgMean,bgStd] = spatialMovAveBG3D(imageLast5,imageSizeX,imageSizeY,imageSizeZ)
% edit of spatialMovAveBG for 3d frames
% EHarry March 2012

% try just running 2d version on each layer
bgMean = NaN(imageSizeX,imageSizeY,imageSizeZ);
bgStd = NaN(imageSizeX,imageSizeY,imageSizeZ);
for z = 1:imageSizeZ
    img = squeeze(imageLast5(:,:,z,:));
    [bgMean_tmp,bgStd_tmp] = spatialMovAveBG(img,imageSizeX,imageSizeY);
    bgMean(:,:,z) = bgMean_tmp;
    bgStd(:,:,z) = bgStd_tmp;
end

