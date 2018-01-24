function c = redGreenColorMap
% REDGREENCOLORMAP creates color map for scores plotting 
%      
%               use this with a followint "colormap(c);" 
%               function call
%
% SYNOPSIS      c = redGreenColorMap
%
% INPUT            :         
% 
% OUTPUT        c  : rgp channels
%                           
% DEPENDENCES   prAlpha uses {                               
%                                       }
%
%               prAlpha is used by { 
%                                           }
%
% Matthias Machacek 07/07/04

r_t = 60;
le  = 0.0;
overlap = 90;

% Blue channel
c(1:256,3)=0.1;

% Red channel
x=(1:128-r_t+overlap);
y=-(1-le)/(128-r_t+overlap) .*x + 1;

% Constant red value
c(1:r_t,1)=1;
% Linear decreasing value
c(r_t+1:128+overlap,1)=y;
% No red value
c(128+overlap+1:256,1)=0;

% Green channel
% No green value
c(1:128-overlap-1,2)=0;

% Linear increasing value
c(128-overlap:256-r_t-1,2)=fliplr(y);
% Constant green value
c(256-r_t:256,2)=1;

c = flipdim(c,1);

