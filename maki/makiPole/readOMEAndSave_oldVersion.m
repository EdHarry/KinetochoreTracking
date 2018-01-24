function [movie,metaData,savePath] = readOMEAndSave_oldVersion( fileAndPath )
% EHarry Feb 2012

[ movie,metaData ] = readOME(fileAndPath);

[path,name] = fileparts(fileAndPath);

savePath = fullfile(path,[name,'.mat']);

saveMovieAndMetaData( movie,metaData,savePath );
%save(savePath,'movie','metaData');


end

