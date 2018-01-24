function [arPARAMK,maPARAMK,TOPOK,bic,errFlag] = armaxCoefKalmanMultiNode(...
       TRAJ,TOPO0,nExoIn,arPARAM0,maPARAM0,tryCONN)
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARMAXCOEFKALMAN fits an ARMAX model to a multi-variate time series      %
%                 which may depend on multiple input time series.         %
%%%%% INPUT DESCRIPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mandatory:
%       TOPO0     : Initial guess for 3D concatenation of 2D adjacency matrices 
%                   for each time lag for connections between ARMA processes        
%   
%
%                   Structure:
%                        *Index in 1st dimension gives node connection
%                        originates from
%                        *Index in 2nd dimension gives node connection goes
%                        to
%                        *Index in 3rd dimension gives time lag for each
%                        connection parameter
%    
%       TRAJ      : 2D Matrix of time series to be fit.
%
%                      Structure:
%                            *Index in 1st dimension gives time point
%                            *Index in 2nd dimension gives node number
%-------------------------------------------------------------------------
%  Optional:
%       nExoIn    : Number of time series in TRAJ which are inputs to the 
%                   network. These must be the first columns in TRAJ.                    
%                   Default: 0
%
%       arPARAM0  : 2D Matrix of autoregressive coefficient initial guesses.
%       maPARAM0  : 2D Matrix of moving average coefficient initial guesses.
%
%                    ARMA parameter matrix structure:
%                         *Index in 1st dimension gives time lag for parameter
%                         *Index in 2nd dimension gives node number
%                   Default: If not included, there are assumed to be no
%                          ARMA parameters in the model.
%
%       tryCONN   : 2D Binary adjacency matrix which determines which connections
%                   present in TOPO0 will be used for the current fit.
%                   Uses the same format as TOPO.
%%%%% OUTPUT DESCRIPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       arPARAMK  : MLE estimation of 2D Matrix of autoregressive coefficients.
%
%       maPARAMK  : MLE estimation of 2D Matrix of moving average coefficients.
%
%       TOPO0     : MLE estimation of 3D concatenation of 2D adjacency matrices 
%                   for each time lag for connections between ARMA processes.
%                   Follows same format as TOPO0.
%       bic       : Bayesian information criterion value for resulting model.
%
%       errFlag   : 0 if execution successful, 1 if error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hunter Elliott %%%   Last updated: January, 2008 %%%%%%%%%%%%%%%%%%%%%
%%%% Adapted from code written by Khuloud Jaqaman  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Output Init.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arPARAMK   = [];
maPARAMK   = [];
TOPOK      = [];
bic        = [];
errFlag    =  0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Input Checks                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% turn off warning:divide by zero
warningState = warning;
warning('off','MATLAB:divideByZero');

%check whether all mandatory input variables were supplied
if nargin < 3
    disp('--armaxCoefKalmanMultiNode: The function requires at least 3 input variables!');
    disp('--Require initial guess topology, time series and # of external inputs.');
    errFlag  = 1;
    return
end

%check for optional inputs and set to default
if nargin < 6
    tryCONN = [];
end
if nargin < 5
    arPARAM0 = [];
    maPARAM0 = [];
end


%get maximum order of X params 
%and # of total elements in network
[nTot1,nTot2,xOrder] = size(TOPO0);
%adjust xOrder for inclusion of current time point
xOrder = xOrder - 1;

if nTot1 ~= nTot2
    disp('TOPO must be square!');
    errFlag = 1;
    return
else
    nTot = nTot1;
end

%Check vaue of nExoIn
if (nExoIn < 0) || (nExoIn >= nTot)
    disp('--armaxCoefKalmanMultiNode: Illegal number of external inputs!');
    errFlag = 1;
    return
else
    %Define # of nodes as #total elements minus external inputs
    nNodes = nTot - nExoIn;
end

%determine if ARMA parameters present
if isempty(arPARAM0) || isempty(maPARAM0)
    disp('--armaxCoefKalmanMultiNode: WARNING No ARMA parameters included!');
    disp('--If node(s) have AR but not MA or vice versa, use NaN.');    
    arma = false;
else
    arma = true;
end

%get highest order of autoregressive coef.
[arOrder,nARp] = size(arPARAM0);

%get highest order of moving average coef.
[maOrder,nMAp] = size(maPARAM0);

%Verify that ARMA param matrices have same # of nodes
if nMAp ~= nARp
    disp('--armaxCoefKalmanMultiNode: Node # disagreement in ARMA parameters!');
    disp('--If a node is to have AR but not MA or vice versa, use NaN.');
    errFlag = 1;
    return
end

%if ARMA params included, verify there are the correct number of nodes
% MA params are implicitly checked by above inequality
if arma && (nARp ~= nNodes)
    disp('--armaxCoefKalmanMultiNode:Incorrect number of nodes in ARMA params!');
    errFlag = 1;
    return
end

%obtain trajectory length and number of trajectories
[trajLength,nTraj,dpth] = size(TRAJ);

%verify that correct number of trajectories were included
if nTot ~= nTraj
    disp('--armaxCoefKalmanMultiNode: Not enough trajectories for # of nodes and external inputs!');
    errFlag = 1;
    return
end

% verify that observational error was included
if dpth ~= 2
    disp('--armaxCoefKalmanMultiNode: Trajectories must include observational error!');
    errFlag = 1;
    return
end

%check if a connectivity was included and if not use all connections
if isempty(tryCONN)
    disp('--armaxCoefKalmanMultiNode: WARNING ! No connectivity detected, using ALL POSSIBLE connections!');
    tryCONN = ones(nNodes,nNodes);
end

%verify that connectivity tryCONN is of proper size
[nChk,nChk2] = size(tryCONN);

if (nChk ~= nNodes) || (nChk2 ~= nNodes)
    disp('--armaxCoefKalmanMultiNode: Input connectivity of incorrect size!');
    errFlag = 1;
    return
end

%if ARMA coef. are present, convert them to partial coef. for minimization
if arma
    for j = 1:nNodes
        arPARAMP0(:,j) = inverseLevDurbExpoAR(arPARAM0(:,j)');
        maPARAMP0(:,j) = inverseLevDurbExpoMA(maPARAM0(:,j)');
    end    
else
    arPARAMP0 = [];
    maPARAMP0 = [];
end

%pad connectivity with connections corresponding to external inputs

%find non-zero external connections in TOPO
tempTOPO = sum(abs(TOPO0),3);
tempTOPO(nExoIn+1:end,nExoIn+1:end) = zeros(nNodes,nNodes);
[exFrom,exTo] = find(tempTOPO);

%Put connections from tryCONN and external inputs into padCONN
padCONN = zeros(nTot,nTot);
if nNodes > 1
    padCONN(nExoIn+1:end,nExoIn+1:end) = tryCONN;
end
for j = 1:nExoIn
    padCONN(exFrom(j),exTo(j)) = 1;
end

%determine number and indices of all connections
[ConnFrom, ConnTo] = find(padCONN);
nConn = length(ConnFrom);

%get an initial estimate of the white noise variance
%excluding external inputs
tmp = [];
for j = (nExoIn+1):nTot
    tmp = cat(1, tmp, TRAJ(:,j,1)); 
end

wnVariance0 = nanvar(tmp); %variance from "previous" iteration
wnVariance = 0.8*wnVariance0; %variance from "current" iteration


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Maximum likelihood estimation of model  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%since in Jones' algorithm the variance of the observational error must be
%divided by the unknown variance of the white noise, the minimization is
%done iteratively until the white noise variance does not change
%significantly between iterations.

%assign initial guess of parameters to paramT for tomlab input

paramT = vectorFromParams(TOPO0,nExoIn,arPARAMP0,maPARAMP0,padCONN,0);

%while the variance changes by more than 5% from one iteration to the next,
%keep regressing
while abs(wnVariance-wnVariance0)/wnVariance0 > 0.05

    %update wnVariance0
    wnVariance0 = wnVariance;

    %divide the observational error by the standard deviation of the white
    %noise in all trajectories
    TRAJ2 = TRAJ;
    wnStd = sqrt(wnVariance);
    for i=1:nTraj
        TRAJ2(:,i,2) = TRAJ2(:,i,2) ./ wnStd;
    end

    %get initial guess of parameters
    param0 = paramT;
    
   % startup tomlab if necessary
    success = startupTomlab;
    if ~success
        error('Tomlab could not be launched!')
    end
    
    prob = conAssign('neg2LnLikelihoodMultiNodeMEXcall',[],[],[],...
        [-10*ones(1,(arOrder+maOrder)*nNodes) -2*ones(1,(xOrder+1)*nConn)],...
        [10*ones(1,(arOrder+maOrder)*nNodes) 2*ones(1,(xOrder+1)*nConn)],...
        'locMinNegLikMEXMN',param0);
    
    prob.PriLevOpt = 0;
    prob.arOrder = arOrder;
    prob.maOrder = maOrder;
    prob.xOrder = xOrder;
    prob.TRAJ = TRAJ2;
    prob.tryCONN = padCONN;
    prob.nExoIn = nExoIn;
    prob.nNodes = nNodes;
    prob.connFrom = ConnFrom;
    prob.connTo = ConnTo;
    
    %minimize -2ln(likelihood) using Tomlab's ucSolve
    % -- 1/25/07 jonas: changed printLevel to 0
    result = tomRun('ucSolve',prob,0,2);



    %proceed if minimization was successful
    if result.ExitFlag == 0
        params = result.x_k';
        proceed = 1;
    else
        proceed = 0;
    end

    %if minimization was successful
    if proceed

        %assign parameter values to paramT to use them as starting guess in
        %next iteration
        paramT = params;
        [arPARAMP,maPARAMP,TOPOK] = paramsFromVector(params',nNodes,arOrder,...
            maOrder,xOrder,padCONN);

        if ~isempty(arPARAMP)
            for j = 1:nNodes
                [arPARAMK(:,j),errFlag] = levinsonDurbinExpoAR(arPARAMP(:,j)');
            end
        else
            arParamK = [];
        end
        if ~isempty(maPARAMP)
            for j = 1:nNodes
                [maPARAMK(:,j),errFlag] = levinsonDurbinExpoMA(maPARAMP(:,j)');
            end
        else
            maParamK = [];
        end


        %obtain likelihood, white noise sequence and white noise variance
        
        
        for n = 1:nNodes

            %go over all trajectories to get innovations and their
            %variances
            %get the innovations, their variances and process white noise
            %using Kalman prediction and filtering
            
            prob.currNode = n;
            % TMP CHANGE !!!!
            [innov,innovVars,wnVec,dum1,dum2,errFlag2(n)] = armaxKalmanInnovMultiNode_SingleNode(TRAJ,n,arPARAMK(:,n)',maPARAMK(:,n)',TOPOK,[],[]);
            wnVarianceSamp(n) = nanmean(innov .^2 ./ innovVars);
            
            [H(n),pVPort(n),errFlag2(n)] = portmanteau(wnVec,10,0.01);
            
            if H(n) == 1
                disp('--armaxCoefKalmanMultiNode: Residuals did not pass portmanteau test!')
                errFlag = 1;
            end
            
            
        end %1:nNodes
        
        if any(errFlag2)
            errFlag = 1;
            disp('--armaxCoefKalmanMultiNode: Error during portmanteau testing!')
        end
        
        neg2LnLikelihoodV = neg2LnLikelihoodMultiNodeMEX(params,prob);

        %calculate mean white noise variance of all trajectories (Eq. 3.14)
        wnVariance = mean(wnVarianceSamp);

        %get number of parameters estimated: arOrder AR coefficients, maOrder MA
        %coefficients, xOrder+1``` X coefficients, and white noise variance
        numParam = nNodes*(arOrder + maOrder) + (xOrder+1)*nConn + nNodes;

        %evaluate the Bayesian Information Criterion
        bic = neg2LnLikelihoodV + log(trajLength)*numParam;

    else %if minimization was not successful

        errFlag = 1;
        return

    end %(if proceed)

end %(while abs(wnVariance-wnVariance0)/wnVariance0 > 0.05)

%%%%%%%%%
%  FIN! %
%%%%%%%%%