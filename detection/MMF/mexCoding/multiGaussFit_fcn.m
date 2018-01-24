function [F,J] = multiGaussFit_fcn(x,image,clusterPixels,psfSigma);

[F,J] = fitNGaussians3D_mexCode_mex(x,image,clusterPixels,psfSigma);

F = norm(F).^2;

J = 2.*sum(J)';

end
