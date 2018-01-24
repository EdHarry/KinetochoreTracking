function dataStruct = makiConvertPoleCoordsToInitCoords( dataStruct )
%MAKICONVERYPOLECOORDSTOINITCOORDS Converts the coords in the
%.poleTransform field to an .initCoord field
% EHarry Jan 2012

if ~isfield(dataStruct,'poleReferenceFrame');
    dataStruct.initCoord = [];
    return
end

if isempty(dataStruct.poleReferenceFrame)
    dataStruct.initCoord = [];
    return
end

if ~isfield(dataStruct,'initCoord');
    dataStruct.initCoord = [];
    return
end

if isempty(dataStruct.initCoord)
    dataStruct.initCoord = [];
    return
end


initCoord = dataStruct.initCoord;

pixelSize = [dataStruct.dataProperties.PIXELSIZE_XY,dataStruct.dataProperties.PIXELSIZE_XY,dataStruct.dataProperties.PIXELSIZE_Z];

%new2OriginalIdxMap(1:length(initCoord)) = struct('map',[]);

for i = 1:length(initCoord)
    poleCoords = dataStruct.poleReferenceFrame(i).poleCoords_cartisian; % get cartisian coords
    
    goodCoordIdx = ~isnan(poleCoords(:,1)); % index to goodCoords (not poles or frames with no poles)
    %new2OriginalIdxMap(iTime).map = find(goodCoordIdx);
    
    allCoord = poleCoords(~isnan(poleCoords(:,1)),:);
    allCoordPix = poleCoords(~isnan(poleCoords(:,1)),1:3) ./ repmat(pixelSize,sum(~isnan(poleCoords(:,1))),1);
    allCoordPix(:,4:6) = deal(0.25);
    amp = initCoord(i).amp;
    amp = amp(~isnan(poleCoords(:,1)),:);
    nSpots = sum(~isnan(poleCoords(:,1)));

    initCoord(i).allCoord = allCoord;
    initCoord(i).allCoordPix = allCoordPix;
    initCoord(i).amp = amp;
    initCoord(i).nSpots = nSpots;
    initCoord(i).map = find(goodCoordIdx);
end

dataStruct.initCoord = initCoord;

end

