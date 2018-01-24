function [steps,nPower]=plotPowerSpectrum(signal,s,freq,draw)
% plotPowerSpectrum calculates and plots a normalized power spectrum 
%
% SYNOPSIS    [steps,nPower]=plotPowerSpectrum(signal,s,freq,draw)
%
% INPUT       signal  : signal vector
%             s       : frames for low-pass filtering in the time domain
%             freq    : frames for low-pass filtering in the frequency domain
%             draw    : [ 0 | 1 ] Turns 'off' or 'on', respectively, the 
%                       plotting of the spectrum
%
% OUTPUT      steps   : normalized frequencies (between 0 and 1)
%             nPower  : normalized power spectrum (only the positive part)
%
% DEPENDENCES plotPowerSpectrum uses { }
%             plotPowerSpectrum is used by {}
%
% Aaron Ponti, 2002

% Check input parameters
if nargin==2
    draw=1;
end
if s<1
    s=1;
end
if ~mod(s,2)
    error('''s'' must be odd');
end

% Low pass filter signal, if needed
if s>1
    s=[-fix(s/2):fix(s/2)]; % Low-pass filtering with gauss
    s=gauss1d(s,1);
    s=s/sum(s);
    signal=conv(signal,s);
%     s=ones(1,s)./s;           % Convolution with a [1/s 1/s 1/s ... ]s vector
%     signal=conv(signal,s);
end

% Calculate frequency step f
l=length(signal)/2;
posP=fix(l)+1;
f=1/l;

% Number of steps
if l==fix(l)
    steps=[0:f:l*f-f];
else
    steps=[0:f:l*f];
end

% Calculate power spectrum
power=fftshift(abs((fft(signal)).^2));

% Crop the positive part of the spectrum
nPower=power(posP:end);

% Normalize
nPower=nPower/max(nPower);

% Filter spectrum if needed
if freq>1
    frq=[-fix(freq/2):fix(freq/2)]; % Low-pass filtering with gauss
    frq=gauss1d(frq,1);
    frq=frq./sum(frq);
    nPower=conv(nPower,frq);
    nPower=nPower(fix(freq/2)+1:end-fix(freq/2));
end

% Draw power spectrum, if needed 
if draw==1
    
    % Display signal
    figure;
    h=plot(signal,'k-');
    title('Signal');
    
    % Display power spectrum
    figure;
    h=plot(steps,nPower,'k-');
    set(h,'LineWidth',2);
    title('Normalized power spectrum');
    xlabel('Normalized frequency');
    ylabel('Normalized power');
    hold off;
        
end
