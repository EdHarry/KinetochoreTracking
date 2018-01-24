function cost = fitTracking( job, usePoleTracks )
% EHarry October 2011
if nargin < 2 || isempty(usePoleTracks)
    usePoleTracks = 1;
else
    usePoleTracks = 0;
end
%
% count=0;
%
% for i = 1:length(job)
%     if usePoleTracks && isfield(job(i).dataStruct,'tracks_pole') && ~isempty(job(i).dataStruct.tracks_pole)
%         count = count + 1;
%         tracksM(count).tracks = job(i).dataStruct.tracks_pole;
%     elseif isfield(job(i).dataStruct,'tracks') && ~isempty(job(i).dataStruct.tracks)
%         count = count + 1;
%         tracksM(count).tracks = job(i).dataStruct.tracks;
%     end
% end
%
% if exist('tracksM','var');
%     cost = fitF2FDMultiple( tracksM );
% else
%     cost = 1e100; % if there was no tracking then introduce a massive penilty
% end

dataStruct = [];

for i = 1:length(job)
    
    dataStruct = [dataStruct,job(i).dataStruct];
    
end

if usePoleTracks
    type = 'pole';
else
    type = 'plate';
end

cost = fitF2FDMultiple( dataStruct,type,1 );

end

