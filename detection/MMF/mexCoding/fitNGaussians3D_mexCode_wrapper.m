function [F,J] = fitNGaussians3D_mexCode_wrapper(x0,image,index,psfSigma,sparseP)

if nargout < 2

	F = fitNGaussians3D(x0,image,index,psfSigma);

else
	[F,J] = fitNGaussians3D(x0,image,index,psfSigma);


J = J.*sparseP;
end
end
