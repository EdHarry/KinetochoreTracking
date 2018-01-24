function [ movie,metaData ] = readOMEMatFile_oldVersion( fileAndPath, timePointsToLoad, channelsToLoad, decon )
% EHarry Feb 2012

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
            eval(['movie(:,:,:,tCount,cCount) = deconMovie_' int2str(ic) '_' int2str(it) ';']);
            clear deconMovie_*
        else
            load(fileAndPath,['movie_' int2str(ic) '_' int2str(it)]);
            if ~exist(['movie_' int2str(ic) '_' int2str(it)],'var')
                error('missing movie data!');
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

