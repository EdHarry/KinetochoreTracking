function [lsmimg,lsm2rgb,fileinfo]=openlsm(param)
%OPENLSM opens a LSM image and convert it into Matlab RGB image
%        according to the vector a conversion lsm2rgb contained
%        in the file lsm2rgb_file. If no file found, then no conversion
%
% SYNOPSIS [lsmimg,lsm2rgb,fileinfo]=openlsm(param);
%	
% INPUT    param.lsm2rgb:  file containig the conversion information
%                          It can be either a '.mat' file saved in the  
%                          image folder and containing a variable called
%                          'lsm2rgb', or the variable itself.
%                          The variable is stored as a vector of conversion,
%                          i.e. lsm2rgb=[2 1 3] -> second channel is red
%                                                  first channel is green
%                                                  third channel is blue
%                          If lsm2rgb_file is not specified, it takes the
%                          default value [1 2 3].
%          param.header:   String containig the header of the open GUI
%          param.filename: String containing the filename of the file to open
%          param.filepath: String containinf the path of the file to open
%                        
% OUTPUT   lsmimg:         lsm image modified or not
%          lsm2rgb:        conversion vector
%          fileinfo:       structure containg .name the name and .path the path
%
% CB, 18-12-01


header='Open LSM file';
lsm2rgb=[1,2,3]; %Default value
filepath=pwd;

if isfield(param,'header')
    header=param.header;
end

if isfield(param,'filepath')
    cd param.filepath;
end

if isfield(param,'filename')
    file=param.filename;
    if file(end-3:end)~='.lsm'
            file=[file,'.lsm'];
    end
    if exist('file')==0
        file=[];
    end
end

if exist('file')==0
    file=[];
end

if isempty(file) 
    [file,filepath]=uigetfile('*.lsm',header);
    cd (filepath);
end

if isfield(param,'lsm2rgb')
    lsm2rgb_file=param.lsm2rgb;
    if isa(lsm2rgb_file,'char')
        if lsm2rgb_file(end-3:end)=='.mat'
            lsm2rgb_file=lsm2rgb_file(1:end-4);
        end
        if exist([lsm2rgb_file,'.mat'])
            load([lsm2rgb_file,'.mat']);
        end;
    end
    if isa(lsm2rgb_file,'double')
        lsm2rgb=lsm2rgb_file;
    end    
end

lsmimg=rtifc(file,1);

if size(lsmimg,3)==3
    lsmimg(:,:,lsm2rgb)=lsmimg;
end;

fileinfo.name=file;
fileinfo.path=filepath;