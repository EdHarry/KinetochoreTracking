function saveMovieAndMetaData_oldVersion( movie,metaData,savePath )
% EHarry Feb 2012

[~,~,~,t,c] = size(movie);

%save(savePath,'metaData');

for ic = 1:c
    for it = 1:t
        eval(['movie_' int2str(ic) '_' int2str(it) '=movie(:,:,:,' int2str(it) ',' int2str(ic) ');']);
    end
end

clear movie ic it

save(savePath);

end

