%function [ gaussList, newAmpBg ] = multiGaussListND_mexCode( maskAmp, maskCoord, X, sigma )
function gaussList = multiGaussListND_mexCode( maskCoord, X, sigma )

gaussList = ones(size(maskCoord,1),size(X,1)+1);


%% GAUSS LIST CALCULATION
for gaussIdx = 1:size(X,1)
    tmp = GaussListND_mexCode(maskCoord,sigma,X(gaussIdx,:));
    gaussList(:,gaussIdx) = tmp(:,1);
end


%% GAUSS FITTING
%I = (G1,G2,G3,....,1) *(a1;a2;a3;...;bg);
%To fit -> a1;a2..;bg = (G1,G2,G3,...,1) \ I

%newAmpBg = gaussList \ maskAmp;


end

