function analysis = poleVSPlateCorr( poleC,plateC,verbose )
%POLEVSPLATECORR corrolations btween plate and pole systems
%   Makes auto and cross correlation plots for pole and plate systems

% INPUT:    poleC - pole coords struct as given by poleCoords.m
%           plateC - (optional) plate coords of the same cells, as given by poleCoords.m
%           verbose - (optional) {0}/1 make plots or not?
% EHarry June 2011

analysis = struct('poleVSPlateSisterCentres',[],'poles',[]);

if nargin < 1 || isempty(poleC)
    error('need at least one input (poleCoords)!');
end

if ~isfield(poleC,'polePoleDistance')
    error('poleC must be a struct returned by poleCoords.m!');
end

if nargin < 2 || isempty(plateC)
    plate = 0; % flag for doing plate comparison
elseif ~isfield(plateC,'sisterList')
    error('plateC must be a struct returned by poleCoords.m');
else
    plate = 1;
end

if nargin < 3 || isempty(verbose)
    verbose = 0;
else
    verbose = 1;
end

%% Analysis

%1: plate vs pole comparisons in sister centre change autoCorr

if plate
    % sister centre autoCorr pole system
    c=1;
    for i = 1:length(poleC)
        for j = 1:2:size(poleC(i).z,1)
            trajPole(c).observations = diff((poleC(i).z(j,:)+poleC(i).z(j+1,:))./2)';
            c=c+1;
        end
    end
    gammaPole = autoCorr(trajPole,20);
    analysis.poleVSPlateSisterCentres.autoCorr.pole = gammaPole;
    
    % sister centre autoCorr plate system
    c=1;
    for i = 1:length(plateC)
        s = plateC(i).sisterList;
        for j = 1:length(s)
            trajPlate(c).observations = diff((s(j).coords1(:,1)+s(j).coords2(:,1))./2);
            c=c+1;
        end
    end
    gammaPlate = autoCorr(trajPlate,20);
    analysis.poleVSPlateSisterCentres.autoCorr.plate = gammaPlate;
end


%2: k-fibre length change corr

c=1;
for i = 1:length(poleC)
    k1= poleC(i).k1(:,:,:,1);
    k2= poleC(i).k2(:,:,:,1);
    k1 = squeeze(sqrt(sum(k1.^2,2)));
    k2 = squeeze(sqrt(sum(k2.^2,2)));
    for j = 1:size(k1,1)
        trajK1(c).observations = diff(k1(j,:))';
        trajK2(c).observations = diff(k2(j,:))';
        c=c+1;
    end
end
gammaK1 = autoCorr(trajK1,20);
gammaK2 = autoCorr(trajK2,20);
gammaKCross = crossCorr(trajK1,trajK2,20);
analysis.poles.kFibreLengthChange.autoCorr.K1 = gammaK1;
analysis.poles.kFibreLengthChange.autoCorr.K2 = gammaK2;
analysis.poles.kFibreLengthChange.crossCorr = gammaKCross;

%3: pole autoCorr in x,y,z
c=1;
for i = 1:length(poleC)
    trajPoleX(c).observations = diff(poleC(i).pole1CoordsRotated(:,1,1));
    trajPoleY(c).observations = diff(poleC(i).pole1CoordsRotated(:,2,1));
    trajPoleZ(c).observations = diff(poleC(i).pole1CoordsRotated(:,3,1));
    c=c+1;
    trajPoleX(c).observations = diff(poleC(i).pole2CoordsRotated(:,1,1));
    trajPoleY(c).observations = diff(poleC(i).pole2CoordsRotated(:,2,1));
    trajPoleZ(c).observations = diff(poleC(i).pole2CoordsRotated(:,3,1));
    c=c+1;
end
gammaPoleX = autoCorr(trajPoleX,20);
gammaPoleY = autoCorr(trajPoleY,20);
gammaPoleZ = autoCorr(trajPoleZ,20);
analysis.poles.poleCorr.x.autoCorr = gammaPoleX;
analysis.poles.poleCorr.y.autoCorr = gammaPoleY;
analysis.poles.poleCorr.z.autoCorr = gammaPoleZ;

%4: poleMSD (rel. to plate)
tracks=NaN(2*length(poleC),size(poleC(1).z,2)*8);
c=1;
for i = 1:length(poleC)
    tracks(c,1:8:end) = poleC(i).pole1CoordsRotated(:,1,1)';
    tracks(c,2:8:end) = poleC(i).pole1CoordsRotated(:,2,1)';
    tracks(c,3:8:end) = poleC(i).pole1CoordsRotated(:,3,1)';
    c=c+1;
    tracks(c,1:8:end) = poleC(i).pole2CoordsRotated(:,1,1)';
    tracks(c,2:8:end) = poleC(i).pole2CoordsRotated(:,2,1)';
    tracks(c,3:8:end) = poleC(i).pole2CoordsRotated(:,3,1)';
    c=c+1;
end
msd = getAllTracksMSqD(tracks,(size(tracks,2)./8)-1);
analysis.poles.msd = msd;

%% plotting

if verbose
    figure % make new figure
    inD=inputdlg({'time lag'}); % ask for time lag
    
    if isempty(inD)
        error('no input!');
    end
    
    lag = str2double(inD{1});
    
    if isnan(lag)
        error('input a number!');
    end
    
    p = 1; %subplot index
    if plate
        subplot(2,2,p); % make plot
        hold on
        plot(lag.*(0:20),analysis.poleVSPlateSisterCentres.autoCorr.plate(:,1),'k');
        myErrorbar(lag.*(0:20),analysis.poleVSPlateSisterCentres.autoCorr.plate(:,1),analysis.poleVSPlateSisterCentres.autoCorr.plate(:,2));
        plot(lag.*(0:20),analysis.poleVSPlateSisterCentres.autoCorr.pole(:,1),'g');
        myErrorbar(lag.*(0:20),analysis.poleVSPlateSisterCentres.autoCorr.pole(:,1),analysis.poleVSPlateSisterCentres.autoCorr.pole(:,2));
        xlabel('time lag');
        ylabel('autoCorr');
        text(lag.*10,0.5,'sister centre autoCorr  ----  plate -> black, pole -> green');
        hold off
        p=p+1;
    end
    subplot(2,2,p); % make plot
    hold on
    plot(lag.*(-20:20),analysis.poles.kFibreLengthChange.crossCorr(:,1),'k');
    myErrorbar(lag.*(-20:20),analysis.poles.kFibreLengthChange.crossCorr(:,1),analysis.poles.kFibreLengthChange.crossCorr(:,2));
    plot(lag.*(-20:20),[NaN(20,1);analysis.poles.kFibreLengthChange.autoCorr.K1(:,1)],'g');
    myErrorbar(lag.*(-20:20),[NaN(20,1);analysis.poles.kFibreLengthChange.autoCorr.K1(:,1)],[NaN(20,1);analysis.poles.kFibreLengthChange.autoCorr.K1(:,2)]);
    plot(lag.*(-20:20),[NaN(20,1);analysis.poles.kFibreLengthChange.autoCorr.K2(:,1)],'r');
    myErrorbar(lag.*(-20:20),[NaN(20,1);analysis.poles.kFibreLengthChange.autoCorr.K2(:,1)],[NaN(20,1);analysis.poles.kFibreLengthChange.autoCorr.K2(:,2)]);
    xlabel('time lag');
    ylabel('corr');
    text(0,0.5,'k-fibre auto and crossCorr  ----  k-fibre crossCorr -> black, K1Auto -> green, k2Auto -> red');
    hold off
    p=p+1;
    
    subplot(2,2,p); % make plot
    hold on
    plot(lag.*(0:20),analysis.poles.poleCorr.x.autoCorr(:,1),'k');
    myErrorbar(lag.*(0:20),analysis.poles.poleCorr.x.autoCorr(:,1),analysis.poles.poleCorr.x.autoCorr(:,2));
    plot(lag.*(0:20),analysis.poles.poleCorr.y.autoCorr(:,1),'g');
    myErrorbar(lag.*(0:20),analysis.poles.poleCorr.y.autoCorr(:,1),analysis.poles.poleCorr.y.autoCorr(:,2));
    plot(lag.*(0:20),analysis.poles.poleCorr.z.autoCorr(:,1),'r');
    myErrorbar(lag.*(0:20),analysis.poles.poleCorr.z.autoCorr(:,1),analysis.poles.poleCorr.z.autoCorr(:,2));
    xlabel('time lag');
    ylabel('autoCorr');
    text(lag.*10,0.5,'pole autoCorr in plate system  ----  x -> black, y -> green, z -> red');
    hold off
    p=p+1;
    
    subplot(2,2,p); % make plot
    hold on
    plot(lag.*(1:((size(tracks,2)/8)-1)),analysis.poles.msd(:,1),'k');
    myErrorbar(lag.*(1:((size(tracks,2)/8)-1)),analysis.poles.msd(:,1),analysis.poles.msd(:,2)./sqrt(analysis.poles.msd(:,3)));
    xlabel('time lag');
    ylabel('msd (\mum.^2)');
    text(lag.*10,median(analysis.poles.msd(:,1)),'pole msd');
    hold off
end
end

