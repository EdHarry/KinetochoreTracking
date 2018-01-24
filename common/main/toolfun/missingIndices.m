function m=missingIndices(v,n)
% missingIndices returns the missing indices in a vector of sequential entries
%
% missingIndices may be used after FIND to extract those row or columns
% indices which are missing.
%
% SYNOPSIS m=missingIndices(v,n)
%
%               | 0 0 1 1 |                            total number of rows
%               | 1 0 1 1 |                                    |
% Example: D =  | 0 0 1 0 |    [y x]=find(D); missingIndices(y,5) returns 5
%               | 1 0 0 0 |                   missingIndices(x,4) returns 2
%               | 0 0 0 0 |                                    |
%                                                      total number of columns
%
% INPUT    v : vector of entries to be checked
%          n : maximum number of entries
%
% OUTPUT   m : missing indices 
%
% Aaron Ponti, March 3rd, 2003

% Make sure v is sorted and does not contain repetitions
v=unique(v); % (unique sorts!)

% make vector with 1 where we have a value, 0 where it's missing
fullVector = zeros(n,1);
fullVector(v) = 1;
% find the missing indices
m = find(~fullVector);

% old version
% % Initialize index vector
% m=[];
% 
% % Find missing entries
% c0=0;
% for c1=1:n
%     if isempty(find(v==c1))
%         c0=c0+1;   
%         m(c0)=c1;
%     end
% end
   