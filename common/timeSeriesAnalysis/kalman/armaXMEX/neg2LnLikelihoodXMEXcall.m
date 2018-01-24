function likelihood = neg2LnLikelihoodXMEXcall(paramV,prob)




sum1 = 0;
sum2 = 0;
numMissing = 0;
totalObs = 0;

% SOME OF THIS STUFF SHOULD BE MOVED TO armaxCoefKalmanMEX OR armaxFitKalmanMEX TO PREVENT IT BEING
% EXECUTED REPEATEDLY DURING MINIMIZATION!!! since now it will be
% different from the non-mex anyways, might as well save the computational time...

%Determine the number of input trajectories
nTraj = length(prob.user.trajOut);

% % Check 1st timeseries for observational error
% if ( nnz(prob.user.trajOut(1).observations(:,2)) ~= 0 )
%     
%     prob.user.obsErrPresent = 1;
% else
%     prob.user.obsErrPresent = 0;
% end
% 
% 
% %Check whether observational error status is the same for all
% %timeseries
% 
% if nTraj > 1
% 
%     for j = 2:nTraj
% 
%         if ( nnz(prob.user.trajOut(j).observations(:,2)) == 0 ) == prob.user.obsErrPresent
% 
%             disp('Warning: Not all timeseries contain observational error! Ignoring those which do and estimating variance!');
% 
%             %If any timeseries are missing observational error, ignore it and
%             %estimate its variance
%             prob.user.obsErrPresent = 0;
% 
%         end
% 
%     end
% 
% end

%Determine the x order
prob.user.xOrder = length(paramV) - (prob.user.arOrder + prob.user.maOrder) - 1;

%Adjust for the observational error variance parameter
if ~prob.user.obsErrPresent
    prob.user.xOrder = prob.user.xOrder - 1;
end

for j = 1:nTraj

   
    trajLength = size(prob.user.trajOut(j).observations,1);
    
    if prob.user.xOrder < 0
        prob.user.trajIn(j).observations = zeros(trajLength,2);   
    end        
    
    %Re-arrange trajectory so that it is stored in the correct order in
    %memory for pointer arithmetic in C
    prob.user.TRAJ = cat(3, prob.user.trajIn(j).observations , prob.user.trajOut(j).observations);
    prob.user.TRAJ = permute( prob.user.TRAJ, [1 3 2]);
    
    %For minimization, wnVariance is assumed to be 1;
    prob.user.wnVariance = 1;
    
    probCall = prob.user;

    [sum1Tmp,sum2Tmp, numMissingTmp] = neg2LnLikelihoodXMEX(paramV,probCall);
    
    sum1 = sum1 + sum1Tmp;
    sum2 = sum2 + sum2Tmp;
    numMissing = numMissing + numMissingTmp;
    totalObs = totalObs + trajLength;    

end

likelihood = (sum1 + (totalObs - numMissing)*log(sum2));

fuckshit = 1.234;