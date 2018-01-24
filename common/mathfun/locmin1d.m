function indx=locmin1d(x)
%LOCMIN1D returns a list of local minima in vector x
%
% SYNOPSIS indx=locmin1d(x)
%
% INPUT x : data vector
%
% OUTPUT indx : index list to local minima in x

% STARTED GD 29-Nov-1999

indx = [];
for( i = 2:(length(x)-1))
   if(x(i)<x(i-1) & x(i)<x(i+1))
      indx = [indx,i];
   end;
end;