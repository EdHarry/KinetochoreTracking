function [ movie,metaData ] = readOMEMatFile( fileAndPath, timePointsToLoad, channelsToLoad, decon, crop )
% EHarry Feb 2012
% [ movie,metaData ] = readOMEMatFile( fileAndPath, timePointsToLoad, channelsToLoad, decon, crop )

if nargin < 2
    timePointsToLoad = [];
end

if nargin < 3
    channelsToLoad = [];
end

if nargin < 4 || isempty(decon)
    decon = 0;
end

switch decon
    case {0,1}
    otherwise
        decon = 0;
end

if nargin < 5 || isempty(crop)
    crop = [];
    doCrop = 0;
else
    doCrop = 1;
end

[~,name,ext] = fileparts(fileAndPath);

if ~strcmp(ext,'.mat')
    error(['readOMEMatFile needs to read a ".MAT" file! ',name,ext,' is not a ".MAT" file!']);
end

load(fileAndPath,'movie_1_1','metaData','t','c');

if ~exist('movie_1_1','var')
    error(['file: ',name,ext,' does not contain movie data!']);
end

if ~exist('metaData','var')
    error(['file: ',name,ext,' does not contain a variable named "metaData"!']);
end

[x,y,z] = size(movie_1_1);

clear movie_1_1

if decon
    vars = whos('-file',fileAndPath);
    decon = ismember('deconMovie_1_1',{vars.name});
end

if decon
    load(fileAndPath,'deconMovie_1_1');
    [x,y,z] = size(deconMovie_1_1);
    clear deconMovie_1_1
end

if isempty(timePointsToLoad)
    timePointsToLoad = 1:t;
end

if isempty(channelsToLoad)
    channelsToLoad = 1:c;
end


if doCrop
    isCrop = any(crop(:,1:3),1);
    if any(isCrop)
        croppedSize = diff(crop(:,1:3),1,1)+1;
        frameSize = [x,y,z];
        croppedSize(1,~isCrop) = frameSize(1,~isCrop);
        % if we don't crop along one of the dimensions, we say we take
        % from 1:end instead
        crop(1,~isCrop) = 1;
        crop(2,~isCrop) = frameSize(1,~isCrop);
        crop = crop(:,1:3);
        % update frame sizes
        x = croppedSize(1);
        y = croppedSize(2);
        z = croppedSize(3);
    else
        % don't crop if not necessary
        doCrop = 0;
    end
end



%movie = movie(:,:,:,timePointsToLoad,channelsToLoad);

movie = NaN(x,y,z,length(timePointsToLoad),length(channelsToLoad));
cCount=0;

for ic = channelsToLoad
    cCount = cCount + 1;
    tCount=0;
    for it = timePointsToLoad
        tCount = tCount + 1;
        if decon
            load(fileAndPath,['deconMovie_' int2str(ic) '_' int2str(it)]);
            if ~exist(['deconMovie_' int2str(ic) '_' int2str(it)],'var')
                error('missing movie data!');
            end
            if doCrop
                eval(['deconMovie_' int2str(ic) '_' int2str(it) ' = deconMovie_' int2str(ic) '_' int2str(it) '(' int2str(crop(1,1)) ':' int2str(crop(2,1)) ',' int2str(crop(1,2)) ':' int2str(crop(2,2)) ',' int2str(crop(1,3)) ':' int2str(crop(2,3)) ');']);
            end
            eval(['movie(:,:,:,tCount,cCount) = deconMovie_' int2str(ic) '_' int2str(it) ';']);
            clear deconMovie_*
        else
            load(fileAndPath,['movie_' int2str(ic) '_' int2str(it)]);
            if ~exist(['movie_' int2str(ic) '_' int2str(it)],'var')
                error('missing movie data!');
            end
            if doCrop
                eval(['movie_' int2str(ic) '_' int2str(it) ' = movie_' int2str(ic) '_' int2str(it) '(' int2str(crop(1,1)) ':' int2str(crop(2,1)) ',' int2str(crop(1,2)) ':' int2str(crop(2,2)) ',' int2str(crop(1,3)) ':' int2str(crop(2,3)) ');']);
            end
            eval(['movie(:,:,:,tCount,cCount) = movie_' int2str(ic) '_' int2str(it) ';']);
            clear movie_*
        end
    end
end

if any(isnan(movie(:)))
    error('not all movie data loaded!');
end

end

