function [x0,lb,ub] = mmfInitGuessLowerUpperBounds3D(maximaPosT,maximaAmpT,...
    bgAmpT,psfSigma,clusterPixels,firstFit)
% Edit of mmfInitGuessLowerUpperBounds to work in 3D, used by detectSpots
%   EHarry 2012

%% ORIGINAL HEADER
% % %MMFINITGUESSLOWERUPPERBOUNDS calculates the initial guess and lower and upper bounds for mixture-model fitting
% % %
% % %SYNOPSIS [x0,lb,ub] = mmfInitGuessLowerUpperBounds(maximaPosT,maximaAmpT,...
% % %    bgAmpT,psfSigma,clusterPixels,firstFit)
% % %
% % %INPUT  maximaPosT : Particle positions.
% % %       maximaAmpT : Particle amplitudes.
% % %       bgAmpT     : Background amplitude.
% % %       psfSigma   : Gaussian sigma for approximating the point spread
% % %                    function.
% % %       clusterPixels: List of pixels belonging to cluster of local maxima.
% % %       firstFit   : 1 if this is the first time a group of local maxima is
% % %                    fitted, 0 otherwise.
% % %       PLEASE SEE detectSubResFeatures2D_V2 FOR PROPER CONTEXT.
% % %
% % %OUTPUT x0         : Initial guess.
% % %       lb         : Lower bound.
% % %       ub         : Upper bound.
% % %
% % %REMARKS This function in its current format is not really written for
% % %general use but as a function to be called by detectSubResFeatures2D_V2.
% % %
% % %Khuloud Jaqaman, August 2011

%feature positions
x0 = maximaPosT; %initial guess
%lb = x0 - 2*psfSigma; %lower bound
lb = x0;
lb(:,1:2) = x0(:,1:2) - 2*psfSigma(1); %lower bound
lb(:,3) = x0(:,3) - psfSigma(2);
minPos = min(clusterPixels);
% lb(lb(:,1)<minPos(1),1) = minPos(1);
% lb(lb(:,2)<minPos(2),2) = minPos(2);
lb(lb(:,1)<minPos(1),1) = minPos(1);
lb(lb(:,2)<minPos(2),2) = minPos(2);
lb(lb(:,3)<minPos(3),3) = minPos(3);
if ~firstFit
    lb(end,:) = minPos;
end
%ub = x0 + 2*psfSigma; %upper bound
ub = x0;
ub(:,1:2) = x0(:,1:2) + 2*psfSigma(1); %upper bound
ub(:,3) = x0(:,3) + psfSigma(2);
maxPos = max(clusterPixels);
% ub(ub(:,1)>maxPos(1),1) = maxPos(1);
% ub(ub(:,2)>maxPos(2),2) = maxPos(2);
ub(ub(:,1)>maxPos(1),1) = maxPos(1);
ub(ub(:,2)>maxPos(2),2) = maxPos(2);
ub(ub(:,3)>maxPos(3),3) = maxPos(3);
if ~firstFit
    ub(end,:) = maxPos;
end

%feature amplitudes
x0 = [x0 maximaAmpT];
% lb(:,3) = eps;
% ub(:,3) = 1;
lb(:,4) = eps;
ub(:,4) = 1;

%background intensity
x0 = x0';
x0 = [x0(:); bgAmpT];
lb = lb';
lb = [lb(:); eps];
ub = ub';
ub = [ub(:); 1];
