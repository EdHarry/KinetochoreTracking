function resp = lineFilter1d(data,d,scales,type)
%LINEFILTER1D simulates different line filters on a profile across the line
%
% SYNOPSIS resp = lineFilter1d(data,d,scales,type)
%
% INPUT data : profile across the line (odd number of samples expected; otherwise
%              zero padding performed)
%       d :    distance between two samples
%       scales : vector with the scales on which filtering is applied
%       type : linetype (1) positive line
%                       (2) negative line
%                       (4) wave line
%
% OUTPUT resp : response with length(scales) x length(data)

% preallocation of the response matrix
resp = zeros(length(scales),length(data));

if(~mod(length(data),2))
   data(length(data)+1) = 0;
end;

switch(type),
case 1, weight = [1,0,-1];
case 2, weight = [-1,0,1];
case 4, weight = [-1,1,-1];
otherwise, error('invalid line type entered');
end;

x = -max(scales)*4:d:max(scales)*4;
xhLength = (length(x)-1)/2;

for(i = 1:length(scales))
   if((type == 1) | (type == 2))
      mb = weight(1)*dogauss1d(x-scales(i),scales(i));
      mf = weight(3)*dogauss1d(x+scales(i),scales(i));
      rb = conv(data,mb);
      rf = conv(data,mf);
      rbC= rb(xhLength:length(data)+xhLength-1);
      rfC= rf(xhLength:length(data)+xhLength-1);
      crfC = clipNegative(rfC);
      crbC = clipNegative(rbC);
      resp(i,:)=scales(i)*sqrt(crbC.*crfC);
   else
      mb = weight(1)*dogauss1d(x-2.0*scales(i),scales(i));
      mf = weight(3)*dogauss1d(x+2.0*scales(i),scales(i));
      mc = weight(2)*dogauss1d(x,scales(i));
      rb = conv(data,mb);
      rf = conv(data,mf);
      rc = conv(data,mc);
      rbC= rb(xhLength:length(data)+xhLength-1);
      rfC= rf(xhLength:length(data)+xhLength-1);
      rcC= rc(xhLength:length(data)+xhLength-1);
      crfC = clipNegative(rfC);
      crbC = clipNegative(rbC);
      crcC = clipNegative(rcC);
      % compute the cubic root
      aux = crcC.*crbC.*crfC;
      for(j = 1:length(aux))
         if(aux(j)>0)
            resp(i,j)=scales(i)*exp(log(aux(j))/3);
         else
            resp(i,j) = 0;
         end;  
      end;
   end;
end;
