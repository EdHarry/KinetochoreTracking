function [innovation,innovationVar,wnVector,stateVec,stateCov,errFlag] = ...
    armaxKalmanInnovMultiNode_SingleNode(TRAJ,iNode,arParam,maParam,TOPO,tryCONN,wnVariance)

%ARMAXKALMANINNOVMULTINODE does forward Kalman prediction and filtering of
%a time series with a variable number of inputs using a Multi-Node ARMAX model.
%
%SYNOPSIS [innovation,innovationVar,wnVector,stateVec,stateCov,errFlag] = ...
%    armaxKalmanInnov(TRAJ,iNode,arParam,maParam,TOPO,tryCONN,wnVariance)
%
%INPUT  TRAJ      : (trajLength x nNodes x 2 ) Matrix of model trajectories 
%                   with measurement uncertainties for trajectroy (:,n,1) in (:,n,2).                  
%                   Missing points should be indicated with NaN. 
%       iNode     : Index of node to be filtered - its trajectory should be
%                   TRAJ(:,iNode,1)
%       arParam   : Autoregressive coefficients for node iNode (row vector).
%       maParam   : Moving average coefficients for node iNode(row vector).
%       TOPO      : (nNodes x nNodes x maxXorder) Matrix of coefficients indicating 
%                   dependence on input.
%       tryCONN   : (nNodes x nNodes ) binary matrix indicating 'current'
%                   connectivity to be used for filtering. Optional.
%                   Default: all connections present in TOPO are used.
%       wnVariance: White noise Variance for node iNode. Optional. Default: 1.
%
%OUTPUT (All outputs are specific to node iNode.)
%       innovation   : Vector of differences between predicted and observed data, or innovations.
%       innovationVar: Vector of innovation variances.
%       wnVector     : Estimated white noise in the process.
%       stateVec     : Forward-predicted state vector
%       stateCov     : Forward-predicted covariance matrix
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%REMARKS The algorithm implemented here is an ARMAX generalized version of
%        the algorithm presented in R. H. Jones,"Maximum Likelihood Fitting 
%        of ARMA Models to Time Series with Missing Observations", 
%        Technometrics 22: 389-395 (1980). All equation numbers used are 
%        those in the paper. However, I do not estimate the observational 
%        error variance, but use that obtained from experimental data or
%        simulated trajectories (and is thus time-dependent).
%        
%Khuloud Jaqaman, January 2006
%Hunter Elliott, March 2007 - Modified to allow multiple nodes/inputs.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

innovation = [];
innovationVar = [];
wnVector = [];
stateVec = [];
stateCov = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used when function was called
if nargin < 5
    disp('--armaxKalmanInnovMultiNode: Not enough input arguments!');
    errFlag = 1;
    return
end

%find trajectory length, # of notes
[trajLength,nNodes,depth] = size(TRAJ);

%Use all connections if tryCONN not supplied
if nargin < 6 || isempty(tryCONN)
    tryCONN = ones(nNodes,nNodes);
end

%assign 1 to WN variance if not supplied
if nargin < 7 || isempty(wnVariance)
    wnVariance = 1;
end

%find and verify size of tryCONN
[sizeCONN1,sizeCONN2] = size(tryCONN);

if sizeCONN1 ~= nNodes || sizeCONN2 ~= nNodes
    disp('size of connectivity "tryCONN" must agree with size of trajectory matrix "TRAJ"')
    errFlag = 1;
    return
end

%insure that observational errors are present in trajectories
if depth ~= 2
    disp('--armaxKalmanInnovMultiNode: Trajectory matrix "TRAJ" should include observational error!');
    errFlag = 1;
    return
end

%make sure that iNode is a valid index
if iNode > nNodes || ~(round(iNode) == iNode)
    disp('--armaxKalmanInnovMultiNode: node index "iNode" must be a valid column index in input trajectory "TRAJ"');
    errFlag = 1;
    return
end

%find arOrder, maOrder and maxXOrder
arOrder = size(arParam,2);
maOrder = size(maParam,2);
maxXOrder  = size(TOPO,3) - 1;

%Extract connections which have non-zero x parameters, are present in
%tryCONN and which connect to iNode.
tmp = sum(abs(TOPO),3) .* tryCONN;
connToiNode = find(tmp(:,iNode));
nConnToiNode = length(connToiNode);

%get maxOrder to determine size of matrices and vectors in Eq. 2.15 - 2.17
maxOrder = max(arOrder,maOrder+1);

%add zeros to ends of arParam and maParam to get vectors of length maxOrder
arParamMod = [arParam zeros(1,maxOrder-arOrder)];
maParamMod = [maParam zeros(1,maxOrder-maOrder)];

TRAJ = cat(1,TRAJ,zeros(maxOrder,nNodes,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation of innovations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%construct matrix F (Eqs. 2.15, 2.16)
transitionMat = diag(ones(maxOrder-1,1),1);
transitionMat(end,:) = arParamMod(end:-1:1); 

%construct column vector G (Eqs. 2.15, 2.16) using the recursions in Eq. 2.13
%note that procErrCov(i) is the covariance of the process at time t 
%with the white noise at time t+i-1, normalized by the white noise variance
%(Eqs. 4.2, 4.3 and 4.4)
procErrCov = ones(maxOrder,1);
for i = 2:maxOrder
    dummy = maParamMod(i-1) + arParamMod(1:i-1)*procErrCov(i-1:-1:1);
    procErrCov(i) = dummy;
end

%calculate wnVariance*G*G'
wnVarianceGG = wnVariance*procErrCov*procErrCov';

%construct row vector H (Eq. 2.17)
observationVec = zeros(1,maxOrder);
observationVec(1) = 1;

%construct matrices of dependence on input
inputCoefMat = zeros(maxOrder,maxXOrder+1,nConnToiNode);
for i = 1:nConnToiNode
    inputCoefMat(end,:,i) = TOPO(connToiNode(i),iNode,end:-1:1);
end

%initialize state vector and its covariance matrix
stateVecT_T = zeros(maxOrder,1); %Z(0|0)
[stateCovMatT_T,errFlag] = covKalmanInit(arParam,maParam,procErrCov,...
    arOrder,maOrder,maxOrder); %P(0|0)
stateCovMatT_T = stateCovMatT_T*wnVariance;

%initialize output variables and indxMiss
innovation = NaN*ones(trajLength,1);
innovationVar = NaN*ones(trajLength,1);
wnVector = NaN*ones(trajLength,1);
stateVec = NaN*ones(maxOrder,trajLength);
stateCov = NaN*ones(maxOrder,maxOrder,trajLength);

%go over all points in trajectory
%note that in the iterations t+1 = timePoint, t = timePoint-1

%initalize vector for inputs
netInputVector = zeros(maxOrder,1);

for timePoint = (maxXOrder+1):trajLength

    for currentInputNumber = 1:nConnToiNode
        %get input vectors
        modxOrder = timePoint - 1 + maxOrder -...
            max(1,timePoint-1+maxOrder-maxXOrder);
        tmpInputVector = TRAJ(timePoint-1+maxOrder-modxOrder:timePoint-1+maxOrder,connToiNode(currentInputNumber),1);
        netInputVector = netInputVector + inputCoefMat(:,:,currentInputNumber) * tmpInputVector;
    end
    
    %predict state at time t+1 given state at time t
    stateVecT1_T = transitionMat*stateVecT_T + netInputVector; %Z(t+1|t), Eq. 3.1
    %reset netInputVector
    netInputVector = zeros(maxOrder,1);
        
    %obtain the predicted state's covariance matrix
    stateCovMatT1_T = transitionMat*stateCovMatT_T*transitionMat' ...
        + wnVarianceGG; %P(t+1|t), Eq. 3.2
    
    if ~isempty(find(isnan(TRAJ(timePoint,connToiNode,1)),1)) %if any connected traj missing this obs
        
        %cannot modify state vector and its covariance matrix predicted 
        %from previous timepoint since there is no observation
        stateVecT_T = stateVecT1_T; %Z(t+1|t+1), Eq. 5.1
        stateCovMatT_T = stateCovMatT1_T; %P(t+1|t+1), Eq. 5.2  
        
        %make sure that covariance matrix is symmetric
        stateCovMatT_T = (stateCovMatT_T+stateCovMatT_T')/2;
        stateVec(:,timePoint) = stateVecT_T;
        stateCov(:,:,timePoint) = stateCovMatT_T;
        
    else %if all connected traj have observations

        %get innovation, dy(t+1) (Eq. 3.8)
        innovation(timePoint) = TRAJ(timePoint,iNode,1) - stateVecT1_T(1);

        %and its variance, V(t+1) (Eq. 3.6 & 3.10)
        innovationVar(timePoint) = stateCovMatT1_T(1,1) + TRAJ(timePoint,iNode,2)^2;

        %calculate delta
        delta = stateCovMatT1_T(:,1)/innovationVar(timePoint); %delta(t+1), Eq. 3.5

        %modify state vector prediction using observation
        stateVecT_T = stateVecT1_T + delta*innovation(timePoint); %Z(t+1|t+1), Eq. 3.4

        %update state covariance matrix
        stateCovMatT_T = stateCovMatT1_T - delta*observationVec*stateCovMatT1_T; %P(t+1|t+1), Eq. 3.7
                
        %make sure that matrix is symmetric
        stateCovMatT_T = (stateCovMatT_T+stateCovMatT_T')/2;
        
        %calculate white noise at this time point using 
        %wn(t+1) = x(t+1|t+1) - x(t+1|t) (Eq. 2.10 with j=1)
        wnVector(timePoint) = stateVecT_T(1) - stateVecT1_T(1);
        
        %store values of forward-predicted state vector and its covariance matrix
        stateVec(:,timePoint) = stateVecT_T;
        stateCov(:,:,timePoint) = stateCovMatT_T;

    end %(if isnan(TRAJinOut(timePoint,1)) ... else ...)
    
end %(for timePoint = 1:trajLength)

return;
%%%%% ~~ the end ~~ %%%%%
