function J = fitNGaussians3D_mexCode_J_wrapper(x0,image,index,psfSigma,sparseP)

J = fitNGaussians3D_mexCode_J_mex(x0,image,index,psfSigma);

%J(abs(J)<1e-4)=0;

%J = J + eps;

J = sparse(J).*sparseP;

end
