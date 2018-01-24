function dataStruct = makiGenerateAlignedTracks( dataStruct )
%MAKIGENERATEALIGNEDTRACKS
%   EHarry Dec 2012

dataStruct = makiGenerateTracks_frameAlignment(dataStruct);
[~,tracksMat_1] = convStruct2MatNoMS(dataStruct.tracks);
dataStruct = makiUpdateClass(dataStruct);
dataStruct = makiAlignAllFrames(dataStruct);

dataStruct = makiGenerateTracks_frameAlignment(dataStruct);
dataStruct = makiUpdateClass(dataStruct);
dataStruct = makiAlignAllFrames(dataStruct);
[~,tracksMat_2] = convStruct2MatNoMS(dataStruct.tracks);

maxCount = 10;
count = 0;
while count < maxCount && (~all(size(tracksMat_1)==size(tracksMat_2)) || ~all(ismember(tracksMat_1,tracksMat_2,'rows')))
    count = count + 1;
    tracksMat_1 = tracksMat_2;
    dataStruct = makiGenerateTracks_frameAlignment(dataStruct);
    dataStruct = makiUpdateClass(dataStruct);
    dataStruct = makiAlignAllFrames(dataStruct);
    [~,tracksMat_2] = convStruct2MatNoMS(dataStruct.tracks);
end

dataStruct.updatedClass = [];
dataStruct.frameAlignment = [];

end

