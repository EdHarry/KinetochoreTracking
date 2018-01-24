function [r,alpha,hMap] = lineHough(magn,ori,smplRate,stdAlpha)
%LINEHOUGH computes the hough transformation for a single straight line
%
% SYNOPSIS [r,alpha,hMap] = lineHough(magn,ori,smplRate,stdAlpha)
%
% INPUT  magn : map with the response from a line / edge filter
%        ori  : map with the local orientation of the line / edge structure
%               if not avalibale set []
%        smplRate : (optional) sampling of the parameter space 
%                   (default: [1 pixel, 1/180 * pi rad])
%        stdAlpha : (optional) standard deviation of the local oriention values
%                   (default: smplRate(2) )
%
% OUTPUT r,alpha : estimated line parameters 
%                  with cos(alpha) * x_1 + sin(alpha) * x_2 - r = 0
%                  for any point x on the line
%        hMap    : full map of the Hough transformation.


if(nargin < 3)
   smplRate = [1,1/180*pi];
end;
if(nargin < 4)
   stdAlpha = smplRate(2);
end;

if(~isempty(ori))
   ori = mod(ori+pi,pi);
end;

dim = size(magn);
ctr = [(dim(2)-1)/2+1,(dim(1)-1)/2+1];

vecA = 0:smplRate(2):pi;
vecR = -sqrt(sum(dim.^2))/2:smplRate(1):sqrt(sum(dim.^2))/2;
if(~mod(length(vecR),2))
   vecR(end+1) = vecR(end)+smplRate(1);
end;
llhVecR = (length(vecR)-1)/2;
hMap = zeros(length(vecA),length(vecR));

auxAlpha = [cos(vecA);sin(vecA)];

for(i = 1:dim(1))
   for( j = 1:dim(2))
      if(magn(i,j)>0)
         rDbl = ( [j,i]-ctr )*auxAlpha;
         if(isempty(ori))
            estA = vecA;
         else
            estA = ones(1,length(vecA))*ori(i,j);
         end;
         w = magn(i,j)*exp(-1/2*(vecA - estA).^2/stdAlpha^2);
         for( k=1:length(vecA))
            hMap(k,round(rDbl(k)/smplRate(1))+llhVecR+1) = ...
               hMap(k,round(rDbl(k)/smplRate(1))+llhVecR+1) + w(k);
         end;
      end;
   end;
end;

[aux,auxI] = max(hMap);
[aux,maxRI]= max(aux);
maxAI = auxI(maxRI);

r = vecR(maxRI);
alpha = vecA(maxAI);


