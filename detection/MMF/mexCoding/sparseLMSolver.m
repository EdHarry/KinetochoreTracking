function [solutionT, residualsT, jacMatT] = sparseLMSolver(options, x0, imageC, clusterPixels, psfSigma, pixelSize, pixelIdx)

%numNZ = length(x0)*length(imageC);

%J = sparse(ones(length(imageC), length(x0)));

%[~, J] = fitNGaussians3D_mexCode_mex(x0,imageC,clusterPixels,psfSigma);

%J = fitNGaussians3D_mexCode_J_wrapper(x0,imageC,clusterPixels,psfSigma);

tic
sparseP = fitNGaussians3D_J_initialSparsePattern(x0, pixelSize, pixelIdx, psfSigma);
disp(['sparsePattern calc time = ' num2str(toc)]);

numNZ = nnz(sparseP);

[~, solutionT] = sparselm('fitNGaussians3D_mexCode_F_wrapper', 'fitNGaussians3D_mexCode_J_wrapper', 'anzp',sparseP ,x0, 0, zeros(size(imageC,1),1), numNZ, -1, 1000, options, imageC,clusterPixels,psfSigma,sparseP);

[residualsT, jacMatT] = fitNGaussians3D_mexCode_mex(solutionT,imageC,clusterPixels,psfSigma);

end


