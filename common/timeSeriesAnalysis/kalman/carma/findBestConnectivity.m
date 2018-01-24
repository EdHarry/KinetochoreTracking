function [bestCONN, allBIC, allarPARAMK, allmaPARAMK, allTOPOK, errMat]...
    =findBestConnectivity(allCONN,TRAJ,networkProp,arPARAM0,maPARAM0,TOPO0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% findBestConnectivity tests CARMA models of a variety of connectivities  %
%                 to see which best describes the input time seriees      %
%%%%% INPUT DESCRIPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mandatory:
%
%       allCONN   : 3D Concatenation of 2D binary adjacency matrices
%
%                        *Index in 1st dimension gives node connection
%                        originates from
%                        *Index in 2nd dimension gives node connection goes
%                        to
%                        These 2D matrices are then concatenated along the
%                        3rd dimension
%
%       TRAJ      : 2D Matrix of time series to be fit.
%
%                      Structure:
%                            *Index in 1st dimension gives time point
%                            *Index in 2nd dimension gives node number
%                            Note: Node number will start with 1 at the
%                            first column which is not an input.
%
%     networkProp:  Structure of network properties, with these fields:
%                        .xOrder : Maximum order of connections between
%                                    ARMA processes.
%                        .arOrder: Maximum order of autoregressive
%                                   processes in the network.
%                        .maOrder: Maximum order of moving average
%                                   processes in the network.
%                
%                       *NOTE* If initial guesses are provided for the
%                       connections, AR or MA parameters, then their order
%                       need not be specified.
%
%-------------------------------------------------------------------------
%  Optional:
%
%     networkProp:  *NOTE* These are optional FIELDS, the structure networkProp 
%                   MUST be input!
%
%                        .externalInputs : Number of time series to be
%                                   considered inputs to the network.
%                                   Default: 0
%
%                        .minGuess : Minimum number of initial guesses to
%                                       try for each connectivity.
%                                       Default: 3
%
%                        .maxGuess : Maximum number of initial guesses to
%                                       try for each connectivity.
%                                       Default: 6
%
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
%       arPARAM0  : 2D Matrix of autoregressive coefficient initial guesses.
%       maPARAM0  : 2D Matrix of moving average coefficient initial guesses.
%
%                    ARMA parameter matrix structure:
%                         *Index in 1st dimension gives time lag for parameter
%                         *Index in 2nd dimension gives node number
%                   Default: If not included, there are assumed to be no
%                          ARMA parameters in the model.
%
%%%%% OUTPUT DESCRIPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       bestCONN  : 2D Binary adjacency matrix which recieved the lowest
%                   BIC.
%
%       allBIC    : Vector of BIC values associated with each connectivity.
%
%       arPARAMK  : MLE estimation of 2D Matrix of autoregressive
%                   coefficients corresponding to the connectivity with lowest BIC.
%
%       maPARAMK  : MLE estimation of 2D Matrix of moving average
%                   coefficients corresponding to the connectivity with
%                   lowest BIC.
%
%       TOPOK     : MLE estimation of 3D concatenation of 2D adjacency matrices 
%                   for each time lag for connections between ARMA
%                   processes for the connectivity of lowest BIC.
%                   Follows same format as TOPO0.
%
%       errMat    : 0 if execution successful, 1 if error.
%                   *Index in 1st dimension is the connectivity number.
%                   *Index in 2nd dimension is the initial guess number.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hunter Elliott %%%   Last updated: January, 2008 %%%%%%%%%%%%%%%%%%%%%
%%%% Adapted from code written by Khuloud Jaqaman  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Output Init.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%check number of input arguments and set defaults when necessary
if nargin < 3
    disp('--findBestConnectivity: Not enough input arguments!');
    disp('--Must at least include connectivity(s), trajectories and network properties!');
    return
end
if nargin < 6
    TOPO0 = [];
    if ~isfield(networkProp,'xOrder')
        disp('--findBestConnectivity: if topology guess not included, xOrder must be!');
        return
    end
end
if nargin < 5
    maPARAM0 = [];
    if ~isfield(networkProp,'maOrder')
        disp('--findBestConnectivity: if AR parameter guess not included, maOrder must be!');
        return
    end
end
if nargin < 4
    arPARAM0 = [];
    if ~isfield(networkProp,'arOrder')
        disp('--findBestConnectivity: if AR parameter guess not included, arOrder must be!');
        return
    end
end

[trajLength,nTraj,dpth] = size(TRAJ);

if dpth < 2
    disp('No observational error in trajectories - assuming to be zero!')
    TRAJ = cat(3,TRAJ,zeros(trajLength,nTraj));
end
if dpth > 2
    disp('--findBestConnectivity: weird trajectory size!');
    return
end  

if isfield(networkProp,'externalInputs')
    nExoInputs = length(networkProp.externalInputs);
else
    nExoInputs = 0;
    networkProp.externalInputs = [];
end

if nExoInputs >= nTraj
    disp('--findBestConnectivity: Too many external inputs!');
    return
end

[nCchk,nCchk2,nConnTry] = size(allCONN);

if (nCchk ~= nCchk2) || ( (nCchk + nExoInputs) ~= nTraj)
    disp('--findBestConnectivity: Wrong size in trajectories/external inputs/connectivities!');
    return
end

%determine number of nodes from trajectory number and external input number
nNodes = nTraj - nExoInputs;

%check if guesses provided for ARMA, and check their size
if ~isempty(arPARAM0)
    [networkProp.arOrder,nARchk] = size(arPARAM0);
    arGuess = true;
else
    arGuess = false;
end

if ~isempty(maPARAM0)
    [networkProp.maOrder,nMAchk] = size(maPARAM0);
    maGuess = true;
else
    maGuess = false;
end

if ~isempty(arPARAM0) && ~isempty(maPARAM0)
    %check node # argreement in ARMA params and set number of nodes
    if (nARchk ~= nMAchk) || (nARchk ~= nNodes)
        disp('Node number disagreement in ARMA parameters!');
        return;
    end
end

%if topology guess is present, check it
if ~isempty(TOPO0)
    [nNchk1,nNchk2,networkProp.xOrder] = size(TOPO0);
    networkProp.xOrder = networkProp.xOrder - 1;
    topoGuess = true;
    %make sure that these inputs do not accept any inputs by checking TOPO
    if ~isempty(find(TOPO0(:,1:nExoInputs,:)))
        disp('Exogenous inputs cannot accept input - check topography!');
        return;
    end
    %verify size agreement
    if (nNodes ~= nNchk1-nExoInputs) || (nNodes ~= nNchk2 - nExoInputs)
        display('Node # disagreement in ext. input number & topology!');
        return
    end
else
    topoGuess = false;
end

%retrieve process order values from structure for convenience
arOrder = networkProp.arOrder;
maOrder = networkProp.maOrder;
xOrder = networkProp.xOrder;

%If no initial guess for any parameters, 
%check if the number of repetitions was included
if ~all([arGuess maGuess topoGuess])
    %retrieve minimum number of repetitions or set to default
    if isfield(networkProp,'minGuess')
       minGuess = networkProp.minGuess;
    else
        minGuess = 3;
    end
    %retrieve maximum number of repetitions or set to default
    if isfield(networkProp,'maxGuess')
       minGuess = networkProp.maxGuess;
    else
        maxGuess = 6;
    end
else
    %if all initial guesses are given, don't repeat regression
    minGuess = 1;
    maxGuess = 1;
end

%initialize memory for regressed parameters and BICs
bestBIC = [Inf,NaN];
allarPARAMK = zeros(arOrder,nNodes,nConnTry) .* NaN;
allmaPARAMK = zeros(maOrder,nNodes,nConnTry) .* NaN;
allTOPOK = zeros(nTraj,nTraj,xOrder+1,nConnTry) .* NaN;
allBIC = zeros(nConnTry,1) .* NaN;
errMat = zeros(nConnTry,maxGuess) .* NaN;

%iterate through all connectivities and calculate best fit and BIC
%for each
for currConn = 1:nConnTry
    
    %Reset success count for next connectivity
    currSuccess = 0;
    
    %Keep trying random guesses until either minGuess successful
    %regressions OR until maxGuess iterations successful or not.
    for currTry = 1:maxGuess
        
        %initialize/clear temporary parameter variables
        tmpARK = zeros(arOrder,nNodes,maxGuess) .* NaN;
        tmpMAK = zeros(maOrder,nNodes,maxGuess) .* NaN;
        tmpTOPOK = zeros(nTraj,nTraj,xOrder+1,maxGuess) .* NaN;
        tmpBIC = zeros(maxGuess,1) .* NaN;
        
        %if no initial AR parameter guess, then create one
        if ~arGuess
            %choose random partial AR params for initial guess
            arPARAMP0 = 4*rand(arOrder,nNodes) - 2;
            %convert these to AR parameters for call to armaxCoefKalmanMultiNode
            for j = 1:nNodes
                arPARAM0(:,j) = levinsonDurbinExpoAR(arPARAMP0(:,j)');
            end
        end

        %if no initial MA parameter guess, then create one
        if ~maGuess
            %choose random partial MA params for initial guess
            maPARAMP0 = 4*rand(maOrder,nNodes) - 2;
            %convert these to MA parameters for call to armaxCoefKalmanMultiNode
            for j = 1:nNodes
                maPARAM0(:,j) = levinsonDurbinExpoMA(maPARAMP0(:,j)');
            end
        end

        %if no initial topology guess, create one
        if ~topoGuess
            TOPO0 = 2*rand(nTraj,nTraj,xOrder+1) - 1;                        
        end
        
        tic;
        
        [tmpARK(:,:,currTry),tmpMAK(:,:,currTry),tmpTOPOK(:,:,:,currTry),...
            tmpBIC(currTry),errFlag] = armaxCoefKalmanMultiNode(...
            TRAJ,TOPO0,nExoInputs,arPARAM0,maPARAM0,allCONN(:,:,currConn));
        
        tSingleReg(currConn,currTry) = toc;
        
        %check if an error occured
        if ~isempty(errFlag)
            %in case a subroutine set errFlag to zero...
            if errFlag ~= 0
                errMat(currConn,currTry) = 1;
            else
                currSuccess = currSuccess + 1;
                errMat(currConn,currTry) = 0;
            end
        else
            %if no error flag, increment the success count
            currSuccess = currSuccess + 1;
            errMat(currConn,currTry) = 0;
        end
        
        %reset errorflag
        errFlag = [];        
        save('findBestConnBACKUP.mat');
        
        %If enough successful regressions completed, exit loop
        if currSuccess >= minGuess
            break
        end
    
    end % 1:numGuess
    
    %determine which guess had best fit and store those parameters
    [minTmpBic,indBestTmpBic] = min(tmpBIC);
    
    allarPARAMK(:,:,currConn) = squeeze(tmpARK(:,:,indBestTmpBic));
    allmaPARAMK(:,:,currConn) = squeeze(tmpMAK(:,:,indBestTmpBic));
    allTOPOK(:,:,:,currConn) = squeeze(tmpTOPOK(:,:,:,indBestTmpBic));
    allBIC(currConn) = tmpBIC(indBestTmpBic);   

    disp(currConn)      

end % 1:nConntry

%determine best connectivity and return it
[bestBic,indBestBic] = min(allBIC);
bestCONN = squeeze(allCONN(:,:,indBestBic));

%%%%%%%%%
%  FIN! %
%%%%%%%%%