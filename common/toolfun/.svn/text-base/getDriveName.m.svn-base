function drive = getDriveName(path)
%getDriveName: Get the drive name of the directory.
%
% SYNOPSIS: drive = getDriveName(path)
%    The input 'path' has to be full path. It can be either Unix or PC format.
%    If no valid 'drive' has been found, [] is returned;

drive = [];
if strcmp(path(1),'/') || strcmp(path(1),'~') %HLE - added case for home directory
   %Consider it as Unix directory.
   fileSepInd = findstr('/',path);
   mntInd = findstr('/mnt/',path);
   if ~isempty(mntInd)
      if length(fileSepInd) > 2
         drive = path(1:fileSepInd(3)-1);
      else
         drive = path(1:end);
      end
   else
      if length(fileSepInd) > 1
         drive = path(1:fileSepInd(2)-1);
      else
         drive = path(1:end);
      end
   end
else
   %Consider it as PC directory.
   fileSepInd = findstr('\',path);
   colonInd = findstr(':',path);
   if colonInd(1) < fileSepInd(1)
      drive = path(1:colonInd);
   end
end
