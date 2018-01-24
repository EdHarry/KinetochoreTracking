function [savedFileName, savePath] = secureSave(varargin)
%SECURESAVE saves without overwriting files
%
% SYNOPSIS savedFileName = secureSave(filename, data,...)
%
% INPUT    filename, data... see save.m
%
% OUTPUT   savedFileName : fileName (without pathname!) of the saved file.
%                        If the file already existed, _1 is added to the
%                        end of the name (or _2 etc). If the file didn't
%                        exist, the input name is used and returned.
%          savePath : Path where the file was saved
%
%c: 10/10/02 dT

fname=varargin{1};
[path,body,nr,ext]=getFilenameBody(fname,'_');
if isempty(path)
    path=pwd;
end;
if isempty(ext)
    ext='.mat';
end;
if isempty(nr)
    fname=[path filesep body ext];
else
    nr = str2double(nr);
    fname=[path filesep body '_' num2str(nr) ext];
end
%if filename already exists, add a (increasing) number at end
while(exist(fname,'file'))
    if isempty(nr)
        nr=1;
    else
        nr=nr+1;
    end;
    fname=[path  filesep body '_' num2str(nr) ext];
end;
vars=sprintf('''%s'',', varargin{2:end});
vars=vars(1:end-1);
savecmd=['save(''' fname ''',' vars ');'];
evalin('caller',savecmd);
if ~isempty(nr)
    savedFileName = [body '_' num2str(nr) ext];
else
    savedFileName = [body ext];
end
savePath = path;
