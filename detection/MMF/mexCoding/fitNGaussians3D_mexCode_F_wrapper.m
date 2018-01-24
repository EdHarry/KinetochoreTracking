function F = fitNGaussians3D_mexCode_F_wrapper(x0,image,index,psfSigma,sparseP)

F = fitNGaussians3D_mexCode_F_mex(x0,image,index,psfSigma);

%C = 100;

%cost = max([lb - x0, zeros(length(x0), 1), x0 - ub], [], 2);

%cost = C.*(norm(cost).^2)./length(F);

%F = sqrt(F.^2 + cost);

end
