function labelPowerSpectrum(steps,nPower,minValue,dt,N)
% labelPowerSpectrum plots and labels the normalized power spectrum returned by plotPowerSpectrum
%
% SYNOPSIS    labelPowerSpectrum(steps,nPower,minValue,dt)
%
% INPUT       signal   : vector of normalized positions as returned by plotPowerSpectrum
%             s        : power spectrum vector as returned by plotPowerSpectrum
%             minValue : value [0..1] below which peak are not labeled
%             dt       : sampling interval (in seconds)
%
% DEPENDENCES labelPowerSpectrum uses { }
%             labelPowerSpectrum is used by {}
%
% Aaron Ponti, 2002

% Checki input parameters
if nargin==3
    time=0;
elseif nargin==4
    N=2*length(nPower);
    time=1;
elseif nargin==5
    time=1;
else
    error('Wrong number of input parameters');
end
    
% Plot the power spectrum
figure;
plot(steps,nPower,'k-');
hold on;

% Only peaks with a power larger than minValue have to be labeled
pos=locmax1d(nPower);
pos=pos(find(nPower(pos)>minValue));

% Go through peaks and label them if their power is > that minValue
for i=1:length(pos)
    if time==1
        if pos(i)~=1 % We don't label the DC
            period=N/pos(i)*dt;
            text(steps(pos(i)),nPower(pos(i)),[num2str(round(period)),'s']);
        end 
    else
        if pos~=1 % We don't label the DC
            text(steps(pos(i)),nPower(pos(i)),num2str(pos(i)));
        end
    end
end
hold off;