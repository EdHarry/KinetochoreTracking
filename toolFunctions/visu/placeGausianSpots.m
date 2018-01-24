function [imageCoords,image] = placeGausianSpots( coords , psfSigma, imageSize)
%PLACEGAUSIANSPOTS places gaussian spots on an image
%
%   coords = nx4 list of spots coords, [x,y,z,A]
%   psfSigma = [psfSigma_xy psfSigma_z]
%   imageSize = [imageSize_x imageSize_y imageSize_z]
%
%   EHarry May 2012

% calulate intensities around the coords of spots

% range is the psf extension range around each spot
range = 3;

%% get relevent image coordinates
imageCoords = [];

for i = 1:size(coords,1)
    minC_x = max(floor(coords(i,1) - range.*psfSigma(1)),1);
    maxC_x = min(ceil(coords(i,1) + range.*psfSigma(1)),imageSize(1));
    minC_y = max(floor(coords(i,2) - range.*psfSigma(1)),1);
    maxC_y = min(ceil(coords(i,2) + range.*psfSigma(1)),imageSize(2));
    minC_z = max(floor(coords(i,3) - range.*psfSigma(2)),1);
    maxC_z = min(ceil(coords(i,3) + range.*psfSigma(2)),imageSize(3));
    
    imageCoords = [imageCoords; makeCube(minC_x,maxC_x,minC_y,maxC_y,minC_z,maxC_z)];
end


if nargout > 1 
    %% get intensities
    psfInteg = zeros(size(imageCoords,1),size(coords,1));
    for i = 1:size(coords,1)
        psfInteg(:,i) = coords(i,4).*GaussListND(imageCoords,[psfSigma(1) psfSigma(1) psfSigma(2)],coords(i,1:3));
    end
    
    
    %% make image
    image = zeros(imageSize(1),imageSize(2),imageSize(3));
    
    % add the intensities
    for i = 1:size(imageCoords,1)
        image(imageCoords(i,1),imageCoords(i,2),imageCoords(i,3)) = sum(psfInteg(i,:));
    end
end
%% SUBFUNCTIONS

% get all pixel coords in a cube
    function cubeCoords = makeCube(xMin,xMax,yMin,yMax,zMin,zMax)
        cubeCoords = NaN(length(xMin:xMax).*length(yMin:yMax).*length(zMin:zMax),3);
        count = 0;
        for y = yMin:yMax
            for z = zMin:zMax
                for x = xMin:xMax
                    count = count + 1;
                    cubeCoords(count,:) = [x,y,z];
                end
            end
        end
    end

end

