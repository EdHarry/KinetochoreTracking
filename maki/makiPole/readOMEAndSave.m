function [movie,metaData,savePath] = readOMEAndSave( fileAndPath,decon )
% EHarry Feb 2012

if nargin < 2 || isempty(decon)
    decon = 0;
end

switch decon
    case {0,1}
        % nothing
    otherwise
        decon = 0;
end

[ movie,metaData ] = readOME(fileAndPath);

[path,name] = fileparts(fileAndPath);

deconMovie = [];

if decon
    deconFile = fullfile(path,[name,'_cmle.r3d']);
    if exist(deconFile,'file') == 2
        deconMovie = readDecon(deconFile,metaData);
    end
end

savePath = fullfile(path,[name,'.mat']);

saveMovieAndMetaData( movie,metaData,savePath,deconMovie );
%save(savePath,'movie','metaData');


end

