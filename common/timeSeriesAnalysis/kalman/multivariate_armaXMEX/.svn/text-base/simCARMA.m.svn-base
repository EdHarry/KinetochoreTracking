function [TRAJ,NOISE] = simCARMA(model,trajLength,noiseSigma,TRAJin)
    
    
% Extract parameters from model 

maPARAM = model.maPARAMS;
TOPOLOGY = model.TOPO;
nNodes = size(TOPOLOGY,2);

%%%%%%%%%%%%%%%%
% INPUT CHECKS %
%%%%%%%%%%%%%%%%

%Check that enough inputs were given
if nargin < 1
    disp('Must include at least a model structure!')
    return;
end

%Check if number of timepoints was input.
if nargin < 2
    trajLength = [];
end

%Check if local noise levels were input.
if nargin < 3
    noiseSigma = [];
end

%Check if exogenous input time series were included
if nargin < 4
    TRAJin = [];
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
    

%Use size of TOPOLOGY matrix to determine number of nodes and
%maximum order of connection coefficients/ AR parameters
[arOrder,nNodes,nNodesChk] = size(TOPOLOGY);
arOrder = arOrder - 1;

%Determine length and number of input time series
[trajLengthIn,nInputs] = size(TRAJin);
if isempty(trajLength)
    if ~isempty(TRAJin)
        nTimePts = trajLength;
        if nInputs >= nNodes
            disp('Too many input trajectories!');
            return
        end
    else
        disp('If no input time series is given, the number of time points to be simulated must be specified!');
        return
    end
elseif ~isempty(TRAJin)
    if nTimePts ~= trajLength
        disp('Input trajectory length does not match # of time points!');
        return
    end
end


%Determine number of sets of ARMA coeficients input

[maOrder,maN] = size(maPARAM);

%Check that the number of nodes in the AR and MA coefficients is the same.
%NOTE: Nodes with only AR or only MA should have zeros in place of missing
%coeeficients
if nNodes ~= maN
    disp('Node # disagreement in ARMA parameters! (For nodes with only AR or only MA use zeros)');
    return;
end

%Verify that 'to' and 'from' dimensions of topology matrix are equal in size
if nNodes ~= nNodesChk
    disp('Topology matrix of improper size!');
    return;
end

%Check if local noise sigmas were input, if not set to zero
if isempty(noiseSigma)
    disp('Assuming no local process noise!');
    noiseSigma = zeros(1,nNodes);
else
    %Check that the correct number of white noise sigmas was given
    if length(noiseSigma) ~= (nNodes-nInputs) || ndims(noiseSigma) ~= 2
        disp('Incorrect number of noise sigma given!');
        return;
    end
end

%check ARMA coef. for causality and invertibility (and stationarity...)
% And make sure that inputs don't have AR parameters
for nNd = 1:nNodes
    
    if nNd <= nInputs
        if (find(TOPOLOGY(:,nNd,:)))
            disp('Input nodes cannot have autogregressive parameters or inputs!');
            return
        end
        if (find(maPARAM(:,nNd)))
            disp('Input nodes cannot have moving average parameters!');
            return
        end
        
    else
        rAR = abs(roots([-squeeze(TOPOLOGY(end:-1:2,nNd,nNd))' 1]));
     
        if ~isempty(find(rAR<=1)) && (nNd > nInputs)
            disp('--simArmaX: Causality requires the AR polynomial not to have any roots for z <= 1!');
            return;
        end
    end
     rMA = abs(roots([maPARAM(end:-1:1,nNd)' 1]));
     
    if ~isempty(find(rMA<=1))
        disp('--simArmaX: Invertibility requires the MA polynomial not to have any roots for z <= 1!');
        return;
    end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trajectory generation                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%allocate memory for trajectories
TRAJ = zeros(trajLength,nNodes);
if (nInputs > 0)
    TRAJ(:,1:nInputs) = TRAJin;
end
NOISE = zeros(trajLength,nNodes);
%initialize normal random number generator
randn('state',sum(100*clock));

%get white noise matrix
for n = nInputs+1:nNodes
    NOISE(:,n) = noiseSigma(n-nInputs) .* randn(trajLength,1);
end

%construct trajectory

for i = max([arOrder maOrder])+1:trajLength
    for ito = 1:nNodes
        
        currSum = 0;
        
        for ifrom = 1:nNodes
            currSum = currSum + squeeze( TOPOLOGY(:,ito,ifrom) )' * TRAJ(i:-1:i-arOrder,ifrom);        % AR / X                
        end
        
        if ito > nInputs
            TRAJ(i,ito) = currSum + maPARAM(:,ito)' * NOISE(i-1:-1:i-maOrder,ito) ...             % MA
                    + NOISE(i,ito);                                                       % noise
        end
        
    end
end


%%%%%%%%%%%%%%%%
%%%   FIN!   %%%
%%%%%%%%%%%%%%%%
