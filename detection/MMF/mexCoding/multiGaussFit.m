function [solutionTAll,residualsT,jacMatT,poptAll,ret,info] = multiGaussFit( x0,image,index,psfSigma )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% [ret, popt, info]=sparselm(fname, jacname, 'an', p0, ncnst, x, jnnz, jtjnnz, itmax, opts, ...)


%maximaPosT = rand(1,3);
%maximaPosT = [100 100 10; 50 50 5; 80 70 4];

%maximaAmpT = rand(1,1);
%maximaAmpT = [1; 1 ; 1];

%bgAmpT = 1;
bgAmpT = 0;

psfSigma = [1 1];

x=100;,y=100;,z=20;,[i,j,k] = ind2sub([x y z],1:x*y*z);
clusterPixels = [i',j',k'];
%[x0,lb,ub] = mmfInitGuessLowerUpperBounds3D(maximaPosT,maximaAmpT,bgAmpT,psfSigma,clusterPixels,1);

%image = rand(x*y*z,1);
image = fitNGaussians3D_mexCode_mex([100.3;100.6;10.4;0.8;95.1;95.8;5.2;0.7;80.6;70.1;3.5;0.9;0.3],zeros(x*y*z,1),clusterPixels,[1 1]);
%image = fitNGaussians3D_mexCode_mex([100.3;100.6;10.4;0.8;0.3],zeros(x*y*z,1),clusterPixels,[1 1]);

cands(1).Lmax = [100 100 10];
cands(1).status = 1;
cands(1).amp = 1;

cands(2).Lmax = [95 95 5];
cands(2).status = 1;
cands(2).amp = 1;

cands(3).Lmax = [80 70 3];
cands(3).status = 1;
cands(3).amp = 1;

clusters = findOverlapPSFs3D_new(cands,x,y,z,psfSigma);

for i = 1:length(clusters)

maximaPosT = clusters(i).maximaPos(:,1:3);
maximaAmpT = clusters(i).maximaAmp;
clusterPixels = clusters(i).pixels(:,1:3);
[x0,lb,ub] = mmfInitGuessLowerUpperBounds3D(maximaPosT,maximaAmpT,bgAmpT,psfSigma,clusterPixels,1);
imageC = image(clusters(i).pixels(:,4));


options = optimset('Jacobian','on',...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-8);
tic
[solutionT,~,residualsT,~,~,~,jacMatT] = ...
            lsqnonlin(@fitNGaussians3D_mexCode_mex,x0,lb,ub,options,imageC,...
            clusterPixels,psfSigma);
time1 = toc
%F = fitNGaussians3D_mexCode_F_mex(x0,image,clusterPixels,psfSigma);
%J = fitNGaussians3D_mexCode_J_mex(x0,image,clusterPixels,psfSigma);

options = [1E-03, 1E-6, 1E-6, 1E-8, 1E-06, sparselm_spsolvr('cholmod')];

%numNZ = length(x0) .* length(image);

% = fitNGaussians3D_mexCode_J_wrapper(x0,image,clusterPixels,psfSigma);
%J = fitNGaussians3D_mexCode_J_wrapper(x0,image,clusterPixels,psfSigma);

%numNZ = nnz(J);
numNZ = length(x0)*length(imageC);

%tic
%optionsM = optimset('Jacobian','on',...
%    'MaxFunEvals', 10, ...
%    'MaxIter', 10, ...
%    'Display', 'off', ...
%    'TolX', 1e-6, ...
%    'Tolfun', 1e-6);

%[solutionTemp,~,~,~,~,~,jacMatTemp] = ...
%            lsqnonlin(@fitNGaussians3D_mexCode_mex,x0,lb,ub,optionsM,image,...
%            clusterPixels,psfSigma);

%time2_1 = toc

%jacMatTemp = nonLnJFun(x0, lb, ub,image,clusterPixels,psfSigma);
%numNZ = nnz(jacMatTemp);



%size(J)

%F = fitNGaussians3D_mexCode_F_wrapper(x0, lb, ub, image, clusterPixels, psfSigma)

tic
%[ret, popt, info] = sparselm('fitNGaussians3D_mexCode_F_mex', 'fitNGaussians3D_mexCode_J_wrapper', 'anzp', J,x0, 0, zeros(size(image,1),1), numNZ, -1, 1000, options, image,clusterPixels,psfSigma);
[ret, popt, info] = sparselm('fitNGaussians3D_mexCode_F_mex', 'fitNGaussians3D_mexCode_J_wrapper', 'an' ,x0, 0, zeros(size(imageC,1),1), numNZ, -1, 1000, options, imageC,clusterPixels,psfSigma);
%[ret, popt, info] = sparselm('fitNGaussians3D_mexCode_F_mex', 'fitNGaussians3D_mexCode_J_wrapper', 'anzp', jacMatTemp, solutionTemp, 0, zeros(size(image,1),1), numNZ, -1, 1000, options, image,clusterPixels,psfSigma);
%[ret, popt, info] = sparselm('nonLnFFun', 'nonLnJFun', 'anzp', jacMatTemp, x0, 0, zeros(size(image,1),1), numNZ, -1, 1000, options, lb, ub,image,clusterPixels,psfSigma);
time2_2 = toc

solutionTAll(i).sln = solutionT;
poptAll(i).popt = popt;

end

end

