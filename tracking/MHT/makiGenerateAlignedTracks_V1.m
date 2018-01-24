function dataStruct = makiGenerateAlignedTracks_V1( dataStruct )
%MAKIGENERATEALIGNEDTRACKS
%   EHarry Dec 2012

dataStruct = makiGenerateTracks_frameAlignment(dataStruct);
nTracks = length(dataStruct.tracks);
dataStruct = makiUpdateClass(dataStruct);
dataStruct = makiAlignAllFrames(dataStruct);

dataStruct = makiGenerateTracks_frameAlignment(dataStruct);
dataStruct = makiUpdateClass(dataStruct);
dataStruct = makiAlignAllFrames(dataStruct);
nNewTracks = length(dataStruct.tracks);

maxCount = 10;
count = 0;
while count < maxCount && nNewTracks < nTracks
    count = count + 1;
    nTracks = nNewTracks;
    dataStruct = makiGenerateTracks_frameAlignment(dataStruct);
    dataStruct = makiUpdateClass(dataStruct);
    dataStruct = makiAlignAllFrames(dataStruct);
    nNewTracks = length(dataStruct.tracks);
end

dataStruct.updatedClass = [];
dataStruct.frameAlignment = [];

end

