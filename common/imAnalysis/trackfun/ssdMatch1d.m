function shift = ssdMatch1d( signal1, signal2, threshold, max_iterations);
%SSDMATCH1D estimates displacement of two signals via SSD
%
% Determines the  displacement of signal1 relative to signal2 with the SSD-Method. 
% Therefore a frame is positioned on signal1. This frame is initially positioned 
% on the same place in signal2. The SSD-Method will shift this frame at signal2 
% in order to match its content to the content of the frame on signal1. In order 
% to obtain high precision this is done by several iteration steps of the SSD-Method. 
%
% SYNOPSIS shift = ssdMatch1d( signal1, signal2, threshold, max_iterations);
%
% An example for calling this function to determine the x-shift between two sinus
% signals: ssd( sin( [ 0: 1 : 360]*pi/180),sin( [ 10 : 1 : 10+360]*pi/180), 0.01, 10).
%
% INPUT signal1: signal at time to
%       signal2: signal at time to + dt
%       threshold: resolution for shift estimaiton 
%       max_iterations: maximum number of iterations  
%
% OUTPUT shift: shift estimate
%
% NOTE sample:
%    ssd( sin( [ 0: 1 : 360]*pi/180),sin( [ 10 : 1 : 10+360]*pi/180), 0.01, 10).
%
%                               
% Robert Muhr, 1998

[m,n] = size(signal1);
x = [ 1 : 1 : n ];
center = round(n/2);    % center of the frame
length = n/2;           % length of the frame
frame = center - length/2 : 1 : center + length/2; 
ref_intensity = interp1( x, signal1, frame, 'cubic');   % reference intensity
new_intensity = interp1( x, signal2, frame, 'cubic');   % new intensity

% intensity derivative of the original signal
dx_matrix = 1/2 * [1, 0, -1];
intensity_derivative = conv( dx_matrix, ref_intensity);

[n,m] = size(ref_intensity);
M = intensity_derivative(1,3:1:m-2)';   % removement of unvalid values
shift = -inv(M'*M)*M'*(new_intensity(1,3:1:m-2)' - ref_intensity(1,3:1:m-2)');

% initialisations
epsilon = threshold + 1;
number_of_iterations = 1;
% iterations start here
while (( abs(epsilon) > threshold ) & (number_of_iterations < max_iterations))
  new_intensity = interp1( x, signal2, frame + shift, 'cubic');
  shift_before =  shift;
  shift = -inv(M'*M)*M'*(new_intensity(1,3:1:m-2)' - ref_intensity(1,3:1:m-2)');
  shift = shift + shift_before;
  number_of_iterations = number_of_iterations + 1;
  epsilon = shift - shift_before;
end;

%display('x-shift of signal2 relative to signal1 determined by the SSD method: ');

%figure(1)
%hold on;
%xlabel('angle (degree)');
%ylabel('Intensity');
%plot( x, signal1, '-b' );
%plot( x, signal2, '-b' );