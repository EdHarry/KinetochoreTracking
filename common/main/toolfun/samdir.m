function success = samdir(dir1,dir2)
%samdir This function checks if two directories are the same.
%
% SYNTAX success = samdir(dir1,dir2)
%    Return 1 if yes, 0 if any directory does not exist or they are not .
%    identical. 'dir1' and 'dir2' are two strings referring to
%    two directories. 

if ~isdir(dir1) | ~isdir(dir2)
   success = 0;
   return;
end

oldDir = pwd;

cd(dir1);
mDir1 = pwd;
cd(dir2);
mDir2 = pwd;

if strcmp(mDir1,mDir2)
   success = 1;
else
   success = 0;
end

cd(oldDir);
