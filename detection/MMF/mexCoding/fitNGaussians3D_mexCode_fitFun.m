function [solutionT,residualsT,jacMatT] = fitNGaussians3D_mexCode_fitFun(options,x0, lb, ub, imageC, clusterPixels, psfSigma, pixelSize, pixelIdx)

sparseP = fitNGaussians3D_J_initialSparsePattern(x0, pixelSize, pixelIdx, psfSigma);

[solutionT,~,residualsT,~,~,~,jacMatT] = lsqnonlin(@fitNGaussians3D_mexCode_wrapper,x0,lb,ub,options,imageC,clusterPixels,psfSigma,sparseP);

end
