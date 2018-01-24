function [fishInfoMat,innovstore] = armaxFisherInfoMatrixMultiNode(TRAJ,...
    arPARAM,maPARAM,TOPO,wnVariance)

tmp = sum(abs(TOPO),3);

[connFrom,connTo] = find(tmp);
nConn = length(connFrom);
%tryCONN = tmp ./ tmp;

[arOrder,nNodes] = size(arPARAM);
[maOrder,nNodesChk] = size(maPARAM);
[nTot,nTotChk,xOrder] = size(TOPO);
[trajLength,nTotChk2,dpth] = size(TRAJ);
xOrder = xOrder-1;
nExoIn = nTot - nNodes;


[paramV,paramI] = vectorFromParams(TOPO,nExoIn,arPARAM,maPARAM,[],1);

numParam = length(paramV);

maxOrder = max(arOrder,maOrder+1);

%pad ARMA parameters so each is of length maxOrder
arPARAMmod = cat(1,arPARAM,zeros( maxOrder-arOrder,nNodes));
maPARAMmod = cat(1,maPARAM,zeros( maxOrder-maOrder,nNodes));

%construct the transition matrix, F
%and the big matrix of F 

FBigMat = [];
FDerivBigMat = [];
GBig = [];
GDerivBig = [];
GGprimeBigMat = [];
GtimesGBigMat = [];
BBigMat = zeros(numParam*maxOrder,numParam*(xOrder+1),nTot);
BDerivBigMat = zeros(numParam*maxOrder,numParam*(xOrder+1),nTot);
HtimesHBigMat = [];
stateCovBigMat = [];
stateCovDerivBigMat = [];

H = zeros(1,maxOrder);
H(1) = 1;
HtimesH = H'*H;

obsIndex = 1:maxOrder:numParam*maxOrder;

for j = 1:numParam

    %determine node associated with current parameter and the
    %connections leading to it
    currNode = paramI.iNode(j);
    iConnToiNode = find(connTo == currNode+nExoIn);
    nConnToiNode = length(iConnToiNode);
    
    %Put numParam copies of F along the diagonal of a big matrix
    F = diag(ones(maxOrder-1,1),1);
    F(end,:) = arPARAMmod(end:-1:1,currNode)';
    FBigMat = blkdiag(FBigMat,F);

    %get derivatives of F as appropriate to parameter type
    %and put along the diagonal of a big matrix
    FDeriv = zeros(maxOrder,maxOrder);    
    if strcmp(paramI.type(j,:),'AR')
       FDeriv(end,maxOrder-paramI.paramLag(j)+1) = 1;
    end
    FDerivBigMat = blkdiag(FDerivBigMat,FDeriv);

    G = ones(maxOrder,1);
    for k = 2:maxOrder
        dummy = maPARAMmod(k-1,currNode) ...
            + arPARAMmod(1:k-1,currNode)'*G(k-1:-1:1);
        G(k) = dummy;
    end
    
    GBig = cat(1,GBig,G);
    
    GDeriv = zeros(maxOrder,1);
    
    if strcmp(paramI.type(j,:),'AR')

        %initialize vectors for the derivatives of G
        paramDeriv = zeros(1,maxOrder);
        paramDeriv(paramI.paramLag(j)) = 1;

        for k = 2:maxOrder
            dummy = paramDeriv(1:k-1)*G(k-1:-1:1) ...
                + arPARAMmod(1:k-1,currNode)'*GDeriv(k-1:-1:1);
            GDeriv(k) = dummy;
        end
   
    elseif strcmp(paramI.type(j,:),'MA')
        
        paramDeriv = zeros(1,maxOrder);
        paramDeriv(paramI.paramLag(j)) = 1;
        
        for k = 2:maxOrder
            dummy = paramDeriv(k-1) + arPARAMmod(1:k-1,currNode)'*GDeriv(k-1:-1:1);
            GDeriv(k) = dummy;
        end
    end %if MA or AR

    GDerivBig = cat(1,GDerivBig,GDeriv);
    
    GGprime = G*G' .* wnVariance(currNode);
    GGprimeBigMat = blkdiag(GGprimeBigMat,GGprime);
    
    GtimesGDeriv = (GDeriv*G'+G*GDeriv') .* wnVariance(currNode);
    GtimesGBigMat = blkdiag(GtimesGBigMat,GtimesGDeriv);
    
    %B uses 'raw' indices like TOPO and TRAJ
    B = zeros(maxOrder,xOrder+1,nTot);
    BDeriv = zeros(maxOrder,xOrder+1,nTot);
    for k = 1:nConnToiNode
        B(end,:,connFrom(iConnToiNode(k))) = squeeze(TOPO(connFrom(iConnToiNode(k)),currNode+nExoIn,end:-1:1));
        if strcmp(paramI.type(j,:),'XX') && connFrom(iConnToiNode(k)) == paramI.fNode(j)+nExoIn
            BDeriv(end,xOrder+1-paramI.paramLag(j),connFrom(iConnToiNode(k))) = 1;
        end
    end
    BBigMat( (j-1)*maxOrder+1:j*maxOrder,(j-1)*(xOrder+1)+1:j*(xOrder+1),:) = B;
    BDerivBigMat( (j-1)*maxOrder+1:j*maxOrder,(j-1)*(xOrder+1)+1:j*(xOrder+1),:) = BDeriv;
    
    HtimesHBigMat = blkdiag(HtimesHBigMat,HtimesH);
    
    %get the initial state covariance matrix
    [stateCovMat00,errFlag] = covKalmanInit(arPARAM(:,currNode)',maPARAM(:,currNode)',...
        G,arOrder,maOrder,maxOrder);
    stateCovMat00 = stateCovMat00 .* wnVariance(currNode);
    stateCovBigMat = blkdiag(stateCovBigMat,stateCovMat00);
    
    %get the partial derivatives of the initial state covariance matrix
    [stateCovMatDeriv00,errFlag] = covKalmanInitDeriv(F,FDeriv,...
        G,GDeriv,wnVariance(currNode),xOrder,stateCovMat00,1e-5);
    
    stateCovDerivBigMat = blkdiag(stateCovDerivBigMat,stateCovMatDeriv00);
    
    
end

fishInfoMat = zeros(numParam);
fishInfoMatA = zeros(numParam);
fishInfoMatB = zeros(numParam);
fishInfoMatC = zeros(numParam);
fishInfoMatD = zeros(numParam);

TRAJ = cat(1,TRAJ,zeros(maxOrder,nTot,2));

%assign initial state vector and it derivatives
stateVecT_T = zeros(maxOrder*numParam,1);
stateVecDerivT_T = zeros(maxOrder*numParam,1);

%assign initial state covariance matrix and its derivatives
stateCovMatT_T= stateCovBigMat;
stateCovMatDerivT_T = stateCovDerivBigMat;

Beta = zeros(numParam,1);

%go over all time points in the trajectory
for j=1:trajLength

                
    %calculate contribution of input series


     for k = 1:nTot
        %get modified xOrder to deal with first few points
        modXOrder = j - 1 + maxOrder - max(1,j-1+maxOrder-xOrder);
                
        %get input vector
        inputVector = TRAJ(j-1+maxOrder-modXOrder:j-1+maxOrder,k,1);
        inputVector = [zeros(xOrder+1-length(inputVector),1); inputVector];
        inputVector = repmat(inputVector,numParam,1);

        %get contribution of input to the state vector and its derivative
        inputContr(:,k) = BBigMat(:,:,k)*inputVector;
        inputContrDeriv(:,k) = BDerivBigMat(:,:,k)*inputVector;
    end

    %predict the state vector
    stateVecT1_T = FBigMat*stateVecT_T + sum(inputContr,2);

    %get the derivative of the predicted state vector
    stateVecDerivT1_T = FDerivBigMat*stateVecT_T + ...
        FBigMat*stateVecDerivT_T + sum(inputContrDeriv,2);

    %obtain the predicted state's covariance matrix
    stateCovMatT1_T = FBigMat*stateCovMatT_T*FBigMat' + GGprimeBigMat;

    %get the partial derivatives of the covariance matrix
    stateCovMatDerivT1_T = FDerivBigMat*stateCovMatT_T*FBigMat' + ...
        FBigMat*stateCovMatDerivT_T*FBigMat' + ...
        FBigMat*stateCovMatT_T*FDerivBigMat' + GtimesGBigMat;

    %update covariance matrix given observation at current time point
    if any(isnan(TRAJ(j,:,1))) %if an observation at this time point is missing

        %"update" state vector and its derivative
        stateVecT_T = stateVecT1_T;
        stateVecDerivT_T = stateVecDerivT1_T;

        %"update" state covariance matrix and its derivatives
        stateCovMatT_T = stateCovMatT1_T;
        stateCovMatDerivT_T = stateCovMatDerivT1_T;

    else %if there is an observation

        %predict observable
        observablesP = stateVecT1_T(obsIndex);

        %calculate innovation
        innovations = TRAJ(j,paramI.iNode+nExoIn)' - observablesP;
        innovstore(:,j) = innovations;
        %get the innovation variance
        innovationVars = diag(stateCovMatT1_T(obsIndex,obsIndex)) + (TRAJ(j,paramI.iNode+nExoIn,2) .^2)';

        %get the partial derivatives of the innovation
        innovDerivs = -stateVecDerivT1_T(obsIndex);

        %get the partial derivatives of the innovation variance
        innovVarDeriv = diag(stateCovMatDerivT1_T);
        innovVarDeriv = innovVarDeriv(obsIndex);    
        
        %calculate delta*H and delta

        tmp = repmat(innovationVars',maxOrder,1);
        innovVarsBigMat = repmat(tmp(:),1,maxOrder*numParam);
        
        deltaH = (stateCovMatT1_T*HtimesHBigMat) ./ innovVarsBigMat;
        delta = deltaH*repmat(H',numParam,1);

        %calculate (derivative of delta) * H and (derivative of delta)
        deltaDerivH = stateCovMatDerivT1_T*HtimesHBigMat ./ innovVarsBigMat ...
            - (stateCovMatT1_T*HtimesHBigMat*stateCovMatDerivT1_T*...
            HtimesHBigMat) ./ (innovVarsBigMat .^2);
        deltaDeriv = deltaDerivH*repmat(H',numParam,1);

        %update the state vector
        innovationsBig = repmat(innovations',maxOrder,1);
        stateVecT_T = stateVecT1_T + delta .* innovationsBig(:);

        %update the derivative of the state vector
        innovDerivsBig = repmat(innovDerivs',maxOrder,1);
        stateVecDerivT_T = stateVecDerivT1_T + deltaDeriv .* innovationsBig(:) ...
            + delta .* innovDerivsBig(:);

        %update state covariance matrix and make sure it's symmetric
        stateCovMatT_T = stateCovMatT1_T - deltaH*stateCovMatT1_T;
        stateCovMatT_T = (stateCovMatT_T+stateCovMatT_T')/2;

        %update derivatives of state covariance matrix
        stateCovMatDerivT_T = stateCovMatDerivT1_T - deltaDerivH*...
            stateCovMatT1_T - deltaH*stateCovMatDerivT1_T;
        stateCovMatDerivT_T = (stateCovMatDerivT_T+stateCovMatDerivT_T')/2;

        %update the Fisher information matrix
        Beta = Beta + innovations.^2 ./ innovationVars;
        


        fishInfoMatA = fishInfoMatA + .25*(innovVarDeriv*innovVarDeriv') ./ (innovationVars*innovationVars');
        fishInfoMatB = fishInfoMatB + ((innovations*innovations') .* (innovDerivs*innovDerivs')) ./ (innovationVars*innovationVars')...
                    + .25*((innovations.^2) * (innovations'.^2) .* (innovVarDeriv*innovVarDeriv')) ./ ((innovationVars.^2)*(innovationVars'.^2))...
                    - .5*((innovations*innovations'.^2) .* (innovDerivs*innovVarDeriv')) ./ (innovationVars*innovationVars'.^2)...
                    - .5*((innovations.^2*innovations') .* (innovVarDeriv*innovDerivs')) ./ (innovationVars.^2*innovationVars');    
        fishInfoMatC = fishInfoMatC + ... 
                    + .5*(repmat(innovations',numParam,1) .* (innovVarDeriv*innovDerivs')) ./ (innovationVars*innovationVars')...
                    - .25*(repmat(innovations'.^2,numParam,1) .* (innovVarDeriv*innovVarDeriv')) ./ (innovationVars*innovationVars'.^2);
        fishInfoMatD = fishInfoMatD + ...        
                    + .5*(repmat(innovations,1,numParam) .* (innovDerivs*innovVarDeriv')) ./ (innovationVars*innovationVars')...
                    - .25*(repmat(innovations.^2,1,numParam) .* (innovVarDeriv*innovVarDeriv')) ./ (innovationVars.^2*innovationVars');
                    
                    

    end %(if isnan(trajOutI(j,1)) ... else ...)

end %(for j=1:trajLength)

fishInfoMatB = fishInfoMatB .* (trajLength)^2 ./ (Beta * Beta');
fishInfoMatC = fishInfoMatC .* trajLength ./ (repmat(Beta',numParam,1));
fishInfoMatD = fishInfoMatD .* trajLength ./ (repmat(Beta,1,numParam));

fishInfoMat = (fishInfoMatA + fishInfoMatB + fishInfoMatC + fishInfoMatD);

%%%%% ~~ the end ~~ %%%%%

%         fishInfoMat = fishInfoMat + (innovVarDeriv*innovVarDeriv') ./ (innovationVars*innovationVars')...
%                     + 4*((innovations*innovations') .* (innovDerivs*innovDerivs')) ./ (innovationVars*innovationVars')...
%                     + ((innovations.^2) * (innovations'.^2) .* (innovVarDeriv*innovVarDeriv')) ./ ((innovationVars.^2)*(innovationVars'.^2))...
%                     + 2*(repmat(innovations',numParam,1) .* (innovVarDeriv*innovDerivs')) ./ (innovationVars*innovationVars')...
%                     - (repmat(innovations'.^2,numParam,1) .* (innovVarDeriv*innovVarDeriv')) ./ (innovationVars*innovationVars'.^2)...
%                     - 2*((innovations*innovations'.^2) .* (innovDerivs*innovVarDeriv')) ./ (innovationVars*innovationVars'.^2)...
%                     + 2*(repmat(innovations,1,numParam) .* (innovDerivs*innovVarDeriv')) ./ (innovationVars*innovationVars')...
%                     - (repmat(innovations.^2,1,numParam) .* (innovVarDeriv*innovVarDeriv')) ./ (innovationVars.^2*innovationVars')...
%                     - 2*((innovations.^2*innovations') .* (innovVarDeriv*innovDerivs')) ./ (innovationVars.^2*innovationVars');
