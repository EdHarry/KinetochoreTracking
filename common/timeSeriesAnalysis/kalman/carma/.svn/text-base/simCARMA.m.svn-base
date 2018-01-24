function [TRAJ,NOISE] = simCARMA(TOPOLOGY,nTimePts,noiseSigma,...
        arPARAM,maPARAM,TRAJin)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simCARMA generates a multivariate Complex Autoregressive Moving Average %
% time series with depedence on input time series.                        %
%%%%% INPUT DESCRIPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       TOPOLOGY  : 3D Concatenation of adjacency matrices at each time lag 
%                   for connections between ARMA processes        
%   
%                   Structure:
%                        *Index in 1st dimension is node connection
%                        originates from
%                        *Index in 2nd dimension is node connection goes
%                        to
%                        *Index in 3rd dimension is time lag for each
%                        connection parameter
%--------------------------------------------------------------------------
%       nTimePts  : Number of time points to be simulated. Not necessary if
%                   input trajector(ies) given.
%--------------------------------------------------------------------------
%       noiseSigma: Standard deviation of white noise at each node. Assumed 
%                   to be normally
%                   distributed with mean zero and given std.
%
%                     Structure:
%                          *Row or column vector of sigma values
%--------------------------------------------------------------------------
%       arPARAM   : Matrix of autoregressive coefficients.
%       maPARAM   : Matrix of moving average coefficients.
%
%                    ARMA parameter matrix structure:
%                         *1st dimension is time lag for parameter
%                         *2nd dimension is node number
%--------------------------------------------------------------------------
%       TRAJin    : 2D Matrix or vector of input time series. If n inputs
%                   are given then the first n nodes in topology matrix are 
%                   assumed to be inputs.
%                    
%                     Structure:
%                           *Longest dimension assumed to be time.
%                           *Remaining dimension is input node
%                           number.
%
%%%%% OUTPUT DESCRIPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       TRAJ      : Matrix of simulated time series
%
%                      Structure:
%                            *1st dimension is time
%                            *2nd dimension is node number
%--------------------------------------------------------------------------
%       NOISE     : Matrix of white noise time series. Structure is same as
%                   TRAJ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hunter Elliott %%%   Last updated: November, 2007  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
% INPUT CHECKS %
%%%%%%%%%%%%%%%%

%Check that enough inputs were given
if nargin < 3
    disp('Must include at least: topology matrix, number of timepoints and local noise sigmas!')
    return;
end

%Check if exogenous input time series were included
if nargin < 6
    TRAJin = [];
    nInputs = 0;
elseif ~isempty(TRAJin)
    %check dimensionality
    if ndims(TRAJin) > 2
        disp('Input trajectory should be 1 or 2D array of timeseries!')
        return
    else
        %Check if 2D
        if ndims(TRAJin) > 1
            %If so then make time the 1st dimension
            if size(TRAJin,1) < size(TRAJin,2)
                TRAJin = TRAJin';
                disp('Assuming Longest dimension in input trajectory is time!')
            end
        end
    end
end


%%% CHECK ARMA COEFICIENTS AND PAD WITH ZEROS %%%

if nargin < 5
    arPARAM = [];
end
if nargin < 4
    maPARAM = [];
end

%if MA are missing but AR are not, fill MA with zeros
if isempty(maPARAM) && ~isempty(arPARAM)
    disp('no moving average (MA) coeficients included.')
    maPARAM = zeros(1,size(arPARAM,2));
end
%If AR are missing but MA is not, fill AR with zeros
if isempty(arPARAM) && ~isempty(maPARAM)
    disp('no autoregressive (AR) coeficients included.')
    arPARAM = zeros(1,size(maPARAM,2));
end

%If no ARMA coef OR input time series then just use size of topology
if nargin < 4 
    arPARAM = zeros(1,size(TOPOLOGY,1));
    maPARAM = arPARAM;
end
%If inputs are included but no ARMA, use topology size adjusted for number
%of inputs
if ~isempty(TRAJin) && (isempty(arPARAM) && isempty(maPARAM))
    arPARAM = zeros(1,size(TOPOLOGY,1)-size(TRAJin,2));
    maPARAM = arPARAM;
end   

%Use size of TOPOLOGY matrix to determine total number of nodes and
%maximum order of connection coefficients
[nTotal,nTotchk,xOrder] = size(TOPOLOGY);

%Determine length and number of input time series
[trajLength,nInputs] = size(TRAJin);
if isempty(nTimePts)
    if ~isempty(TRAJin)
        nTimePts = trajLength;
    else
        disp('If no input time series is given, the number of time points to be simulated must be specified!');
        return
    end
elseif ~isempty(TRAJin)
    if nTimePts ~= trajLength
        disp('Input trajectory length does not match # of time points!')
        return
    end
end


%Determine number of sets of ARMA coeficients input
[arOrder,arN] = size(arPARAM);
[maOrder,maN] = size(maPARAM);
% (connection order of -1 is no connection)
xOrder = xOrder-1;

%Check that the number of nodes in the AR and MA coefficients is the same.
%NOTE: Nodes with only AR or only MA should have zeros in place of missing
%coeeficients
if arN ~= maN
    disp('Node # disagreement in ARMA parameters! (For nodes with only AR or only MA, fill with zeros!)');
    return;
else
    nNodes = arN;
end

% Verify that number of nodes in ARMA coef. is same as Number of
% trajectories and inputs
if nNodes ~= (nTotal - nInputs) | nNodes ~= (nTotchk - nInputs)
    disp('Topology size and number of inputs/nodes disagree!');
    return
end

%Check if local noise sigmas were input, if not set to zero
if isempty(noiseSigma)
    disp('Assuming no local process noise!')
    noiseSigma = zeros(1,nNodes);
else
    %Check that the correct number of white noise sigmas was given
    if length(noiseSigma) ~= nNodes | ndims(noiseSigma) ~= 2
        disp('Incorrect number of noise sigma given!');
        return
    end
end

%check ARMA coef. for causality and invertibility (and stationarity...)
for nNd = 1:nNodes
    
     rAR = abs(roots([-arPARAM(end:-1:1,nNd)' 1]));
     
    if ~isempty(find(rAR<=1))
        disp('--simArmaX: Causality requires the AR polynomial not to have any roots for z <= 1!');
        return;
    end
   
     rMA = abs(roots([maPARAM(end:-1:1,nNd)' 1]));
     
    if ~isempty(find(rMA<=1))
        disp('--simArmaX: Invertibility requires the MA polynomial not to have any roots for z <= 1!');
        return;
    end
  
end

%Check that the nodes corresponding to exogenous inputs do not themselves
%recieve input
if ~isempty(find(TOPOLOGY(:,1:nInputs,:)))
    disp('Inputs cannot accept input - check topology!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trajectory generation                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%allocate memory for trajectories
TRAJ = zeros(nTimePts,nTotal);

%initialize normal random number generator
randn('state',sum(100*clock));

%get white noise matrix
noiseSigma = [ zeros(1,nInputs) noiseSigma];
for n = 1:nTotal
    NOISE(:,n) = noiseSigma(n) .* randn(nTimePts,1);
end

if nInputs > 0
    TRAJ(:,1:nInputs) = TRAJin;  %First trajectories are inputs
end

%construct trajectory

for i = max([arOrder maOrder xOrder ])+1:nTimePts
    for j = nInputs+1:nTotal  %starts at nInputs+1 - inputs aren't ARMA
        TRAJ(i,j) = arPARAM(:,j-nInputs)' * TRAJ(i-1:-1:i-arOrder,j) ...        % AR
            + maPARAM(:,j-nInputs)' * NOISE(i-1:-1:i-maOrder,j) ...             % MA
            + NOISE(i,j);                                                       % noise
            for k = 1:nTotal
                for l = 1:xOrder+1
                    TRAJ(i,j) = TRAJ(i,j) + (TOPOLOGY(k,j,l) * TRAJ(i-l+1,k));  % X
                end
            end
    end
end


%%%%%%%%%%%%%%%%
%%%   FIN!   %%%
%%%%%%%%%%%%%%%%
