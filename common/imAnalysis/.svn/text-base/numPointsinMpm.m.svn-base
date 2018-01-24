function[np]=numPointsinMpm(mpm);
% numPointsinMpm calculates number of points in each frasme of mpm
%
% SYNOPSIS   [np]=numPointsinMpm(mpm);
[sx,sy]=size(mpm);
for i=1:(sy/2)
    ct=1+(i-1)*2;
    xvec=mpm(:,ct);
    np(i)=length(nonzeros(xvec));
end