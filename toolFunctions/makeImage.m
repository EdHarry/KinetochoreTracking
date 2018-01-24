function image = makeImage( imSize, spots, bg, psfSigma )
%MAKEIMAGE makes a fake image of subres 3d gaussians
% imSize -> [x,y,z] size of image in pixels
% spots -> [x,y,z,A;...] matrix of spot position and amplitudes
% bg -> background value
% psfSigma -> [sigmaXY,sigmaZ] -> size of gaussian (in pixels)
%   EHarry Nov 2012

% make 1d list of spot pos
x0 = spots';
x0 = x0(:);
x0(end+1) = bg;

% get complete pixel index list
[x,y,z] = ind2sub(imSize,1:prod(imSize));
pixels = [x' y' z'];

% use fitNGaussians to generate the image
image = fitNGaussians3D_mexCode_mex(x0,zeros(prod(imSize),1),pixels,psfSigma);

% reshape image
image = reshape(image',imSize);

end

