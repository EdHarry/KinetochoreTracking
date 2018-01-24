function J = fitNGaussians3D_J_initialSparsePattern(x0, pixelSize, pixelIdx, psfSigma)
% EHarry, Oct 2012
% x0 -> initial input vector [x1;y1;z1;A1;x2...;bg]
% pixelSize -> size of image [X,Y,Z]
% pixelIdx -> linear (single dimensional) index to the pixels of the image
%              that will be used durin the fitting 
% psfSigma -> gaussian STD [xy,z]   

% range around each spot, this should be identical to the range when clustering spots 
tmpSize = round(9*psfSigma);
tmpSize = tmpSize + (1-mod(tmpSize,2));
psfRange = floor(tmpSize/2);

numPix = length(pixelIdx);
numCol = length(x0);
numK = (numCol-1)/4;

% initialise the pattern to an empty sparse matrix, and make the last
% column all = to 1 (for the background)
J = sparse(numPix, numCol);
J(:,end) = 1;

% inital fitting image, equal to the size of the image plus 2*psfRange to take
% care of spots on the edge 
im = false(pixelSize(1)+2*psfRange(1), pixelSize(2)+2*psfRange(1), pixelSize(3)+2*psfRange(2));

% loop over spots
for k = 1:numK
    x = ceil(x0(((k-1)*4)+1)); % get coords of the spot
    y = ceil(x0(((k-1)*4)+2));
    z = ceil(x0(((k-1)*4)+3));
    
    %disp(['psfRange = ' int2str(psfRange) ', [x,y,z] = ' int2str([x y z])]);
    
    % range of 2 psfRange around each spot (note that ranges are shifted by 1 psfRange because of the extra large image)
    ymin = y; 
    ymax = y+2*psfRange(1);
    xmin = x;
    xmax = x+2*psfRange(1);
    zmin = z;
    zmax = z+2*psfRange(2);
    
    % copy the fitting image and set the range around the spot to true
    im_tmp = im;
    im_tmp(xmin:xmax, ymin:ymax, zmin:zmax) = true;
    
    % crop down the image to its correct size
    im_tmp = im_tmp(psfRange(1)+1:end-psfRange(1),psfRange(1)+1:end-psfRange(1),psfRange(2)+1:end-psfRange(2));
    
    % find the connected cluster in the image (there is only one) and get
    % its linear pixel indexes
    CC = bwconncomp(im_tmp);
    idx = CC.PixelIdxList{1};
    
    % find which rows of the jacobian (pixels) should be used for this spot 
    rowIdx = ismember(pixelIdx, idx);
    
    % make those members = 1 in the jacobian pattern
    J(rowIdx, (((k-1)*4)+1) : (((k-1)*4)+4)) = 1;
end
end
