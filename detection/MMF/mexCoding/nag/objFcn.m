function [mode, f, fjac, user] = objFcn(mode, m, n, ldfj, needfi, x, fjac, nstate, user)%#codegen

% defaults
f = [];
fjac = [];

%% core
%extract background intensity from x0 and remove from vector
bgAmp = x(end);
x = x(1:end-1);

%reshape 4nx1 vector x0 into nx4 matrix
n = (n-1)/4;
x = reshape(x,4,n);
x = x';

%extract PSF center positions and amplitudes
psfPos = x(:,1:3);
psfAmp = x(:,4);

%% index mod
index = user.index;
psfSigma = user.psfSigma;

if needfi > 0
    index = index(needfi,:);
end

%% pixel integration calculations
%find minimum and maximum pixel indices
minIndxX = min(index(:,1));
maxIndxX = max(index(:,1));
minIndxY = min(index(:,2));
maxIndxY = max(index(:,2));
minIndxZ = min(index(:,3));
maxIndxZ = max(index(:,3));

% get xy-indices relative to minimum
% relIndxX = index(:,1) - minIndxX + 1;
% relIndxY = index(:,2) - minIndxY + 1;
% get xy-indices relative to minimum
relIndxX = index(:,1) - minIndxX + 1;
relIndxY = index(:,2) - minIndxY + 1;
relIndxZ = index(:,3) - minIndxZ + 1;

%determine the contribution of each PSF (assuming amplitude 1) to a
%pixel based on its x-coordinate (needed to calculate F & J)
psfIntegX = zeros(maxIndxX-minIndxX+1,n);
for i=1:n
    temp = GaussListND_mexCode((minIndxX:maxIndxX)',...
        psfSigma(1),psfPos(i,1));
    
    temp2 = squeeze(temp);
    %clear temp
    
    psfIntegX(:,i) = temp2(:,1);
end

%determine the contribution of each PSF (assuming amplitude 1) to a
%pixel based on its y-coordinate (needed to calculate F & J)
psfIntegY = zeros(maxIndxY-minIndxY+1,n);
for i=1:n
    temp = GaussListND_mexCode((minIndxY:maxIndxY)',...
        psfSigma(1),psfPos(i,2));
    
    temp2 = squeeze(temp);
    %clear temp
    
    psfIntegY(:,i) = temp2(:,1);
end

%determine the contribution of each PSF (assuming amplitude 1) to a
%pixel based on its z-coordinate (needed to calculate F & J)
psfIntegZ = zeros(maxIndxZ-minIndxZ+1,n);
for i=1:n
    temp = GaussListND_mexCode((minIndxZ:maxIndxZ)',...
        psfSigma(2),psfPos(i,3));
    
    
    temp2 = squeeze(temp);
    %clear temp
    
    psfIntegZ(:,i) = temp2(:,1);
end

%% PSF values for J
if mode > 0
    %calculate the value of each PSF (assuming amplitude 1) at the
    %x-coordinates of the corners of all pixels (needed to calculate J)
    psfValueX = zeros(maxIndxX-minIndxX+2,n);
    for i=1:n
        psfValueX(:,i) = exp(-((minIndxX-0.5:maxIndxX+0.5)'...
            -psfPos(i,1)).^2/2/psfSigma(1)^2);
    end
    
    %calculate the value of each PSF (assuming amplitude 1) at the
    %y-coordinates of the corners of all pixels (needed to calculate J)
    psfValueY = zeros(maxIndxY-minIndxY+2,n);
    for i=1:n
        psfValueY(:,i) = exp(-((minIndxY-0.5:maxIndxY+0.5)'...
            -psfPos(i,2)).^2/2/psfSigma(1)^2);
    end
    
    %calculate the value of each PSF (assuming amplitude 1) at the
    %z-coordinates of the corners of all pixels (needed to calculate J)
    psfValueZ = zeros(maxIndxZ-minIndxZ+2,n);
    for i=1:n
        psfValueZ(:,i) = exp(-((minIndxZ-0.5:maxIndxZ+0.5)'...
            -psfPos(i,3)).^2/2/psfSigma(2)^2);
    end
    
    %% j calculation
    fjac = ones(m,4*n+1); %(last column for background amplitude)
    fjac(:,1:4:4*n) = repmat(psfAmp',m,1).*(psfValueX(relIndxX,:)-...
        psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:).*psfIntegZ(relIndxZ,:); %w.r.t. x
    fjac(:,2:4:4*n) = repmat(psfAmp',m,1).*(psfValueY(relIndxY,:)-...
        psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:).*psfIntegZ(relIndxZ,:); %w.r.t. y
    fjac(:,3:4:4*n) = repmat(psfAmp',m,1).*(psfValueZ(relIndxZ,:)-...
        psfValueZ(relIndxZ+1,:)).*psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); %w.r.t. z
    fjac(:,4:4:4*n) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:).*psfIntegZ(relIndxZ,:); %w.r.t. amp
end


%% f calulation
if mode ~= 1
    if needfi > 0
        f = zeros(m,1);
        f(needfi) = (sum(repmat(psfAmp,1,m).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)'.*psfIntegZ(relIndxZ,:)',1))' + repmat(bgAmp,m,1);
    else
        f = (sum(repmat(psfAmp,1,m).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)'.*psfIntegZ(relIndxZ,:)',1))' + repmat(bgAmp,m,1);
    end
end



end

