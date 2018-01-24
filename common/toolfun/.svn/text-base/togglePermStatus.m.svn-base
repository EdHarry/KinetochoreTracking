function imgOut = togglePermStatus(imgIn)
%TOGGLEPERMSTATUS changes the permutation between 'C' and 'M'
%
%	SYNOPSIS img = togglePermStatus(img)
%
%  INPUT imgIn : image structure with at least the two fields 
%                *.data : 2d image data matrix
%                *.perm : permutation status
%
%  OUTPUT imgOut : the same image structure as imgIn with toggled
%                  permutation status

imgOut = imgIn;
if(~isstruct(imgIn))
   return;
end;
if(~isfield(imgIn,'data') | ~isfield(imgIn,'perm'))
   return;
end;
if(~ndims(imgIn) == 2)
   return;
end;

imgOut.data = permute(imgIn.data,[2,1]);

switch imgIn.perm
case 'M', imgOut.perm = 'C';
otherwise imgOut.perm = 'M';
end;

   
      