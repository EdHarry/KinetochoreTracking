function saveMovieAndMetaData( movie,metaData,savePath,deconMovie )
% EHarry Feb 2012

[~,~,~,t,c] = size(movie);

%save(savePath,'metaData');

for ic = 1:c
    for it = 1:t
        eval(['movie_' int2str(ic) '_' int2str(it) '=movie(:,:,:,' int2str(it) ',' int2str(ic) ');']);
    end
end

if ~isempty(deconMovie)
    [~,~,~,t,c] = size(deconMovie);
    for ic = 1:c
        for it = 1:t
            eval(['deconMovie_' int2str(ic) '_' int2str(it) '=deconMovie(:,:,:,' int2str(it) ',' int2str(ic) ');']);
        end
    end
end

clear movie ic it deconMovie

save(savePath);

end

