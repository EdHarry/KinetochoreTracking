function c=minus(a,b)
% This functions implements the minus operator for variables of class int16
%
% Aaron Ponti, May 7th, 2004

% Min and max values allowed
cmin=-32768;
cmax=32767;

% Perform the summation operation
c=double(a)-double(b);

% Check for possible overflows
if any(c(:)<cmin) | any(c(:)>cmax)
    error('Result out of range.')
else
    c=int16(c);
end
