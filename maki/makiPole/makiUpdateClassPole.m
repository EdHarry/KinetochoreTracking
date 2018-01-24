function dataStruct = makiUpdateClassPole( dataStruct )
%MAKIUPDATECLASSPOLE updates the class and phase of a cell ignoring the
%poles, which otherwise would incorrectly be identified as (outling) kinetochores
% EHarry Jan 2012

dataStruct_temp = dataStruct; % copy the dataStruct for the operations

dataStruct_temp.planeFit = []; % delete the original planeFit

dataStruct_temp = makiConvertPoleCoordsToInitCoords(dataStruct_temp); % convert the poleCoords to initCoords

try
    dataStruct_temp = makiFitPlane(dataStruct_temp,0); % try to fit a new plane
catch
    dataStruct.planeFit_pole = [];
    return
end

% now copy the new phase, unaligned idx and lagging idx

planeFit_new = dataStruct_temp.planeFit;
initCoord_new = dataStruct_temp.initCoord;
planeFit = dataStruct.planeFit;

for i = 1:length(planeFit)
    planeFit(i).phase = planeFit_new(i).phase;
    planeFit(i).unalignedIdx = initCoord_new(i).map(planeFit_new(i).unalignedIdx)';
    planeFit(i).laggingIdx = initCoord_new(i).map(planeFit_new(i).laggingIdx)';
end



%% dependencies
% planeFit_pole depends on dataProperites, initCoords, planeFit, tracks,
% sisterList, updatedClass, frameAlignment, poles and poleReferenceFrame
% dependencies = struct('dataProperties',[],'initCoord',[],'planeFit',[],'tracks',[],'sisterList',[],'updatedClass',[],'frameAlignment',[],'poles',[],'poleReferenceFrame',[]);

dataPropName = dataStruct.dataPropertiesName;
dataPropV = getVersion(dataPropName);
initCoordName = dataStruct.initCoordName;
initCoordV = getVersion(initCoordName);
planeFitName = dataStruct.planeFitName;
planeFitV = getVersion(planeFitName);
tracksName = dataStruct.tracksName;
tracksV = getVersion(tracksName);
sisterListName = dataStruct.sisterListName;
sisterListV = getVersion(sisterListName);
updatedClassName = dataStruct.updatedClassName;
updatedClassV = getVersion(updatedClassName);
frameAlignmentName = dataStruct.frameAlignmentName;
frameAlignmentV = getVersion(frameAlignmentName);
polesName = dataStruct.polesName;
polesV = getVersion(polesName);
poleReferenceFrameName = dataStruct.poleReferenceFrameName;
poleReferenceFrameV = getVersion(poleReferenceFrameName);

dependencies.dataProperties = dataPropV;
dependencies.initCoord = initCoordV;
if ~isempty(dataStruct.planeFit)
    dependencies.planeFit = planeFitV;
end
if ~isempty(dataStruct.tracks)
    dependencies.tracks = tracksV;
end
if ~isempty(dataStruct.sisterList)
    dependencies.sisterList = sisterListV;
end
if ~isempty(dataStruct.updatedClass)
    dependencies.updatedClass = updatedClassV;
end
if ~isempty(dataStruct.frameAlignment)
    dependencies.frameAlignment = frameAlignmentV;
end
dependencies.poles = polesV;
dependencies.poleReferenceFrame = poleReferenceFrameV;

% save into the first planeFit_pole
planeFit(1).dependencies = dependencies;




dataStruct.planeFit_pole = planeFit; % save as a new struct

end

