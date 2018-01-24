function [img,info] = rawread(filename)
% RAWREAD reads raw data image
%
% SYNOPSIS [img,info] = rawread(filename)
%
% INPUT filename: filename string (of .raw file)
%
% OUTPUT img : image
%        info: image info (size, bit depth)

if isempty(findstr(filename,'.'))
	filename=[filename,'.raw'];
end;

[file, message] = fopen(filename,'r');
if file == -1
  a=['file ',filename,' not found.'];
  error(a);
end
info=fread(file,18,'int16');
img=fread(file,[1300,1030],'int16');
img=img';
fclose(file);
