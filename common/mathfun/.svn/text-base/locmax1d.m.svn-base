function indx=locmax1d(x)
%LOCMAX1D returns a list of local maxima in vector x
%
% SYNOPSIS indx=locmax1d(x)
%
% INPUT x : data vector
%
% OUTPUT indx : index list to local maxima in x

% STARTED GD 29-Nov-1999

indx = [];
for( i = 2:(length(x)-1))
   if(x(i)>x(i-1) & x(i)>x(i+1))
      indx = [indx,i];
   end;
end;