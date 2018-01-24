function dataStruct = makiMakeAlignedSisterList( dataStruct )
%MAKIMAKEALIGNEDSISTERLIST Makes a sisterList_pole based only on alignedCoords from a pole based system and adds it to the dataStruct
% EHarry Jan 2012

try
   sisterList_new = makiConstructAlignedSisters_pole(dataStruct); % make the aligned sisterList 
catch
    dataStruct.sisterList_pole_aligned = [];
    return
end

if isempty(sisterList_new)
    dataStruct.sisterList_pole_aligned = [];
    return
end

[sisterList_new.coords1] = sisterList_new.coords1Aligned; % copy the aligned coords to the regualar coords
[sisterList_new.coords2] = sisterList_new.coords2Aligned;

[sisterList_new.distances] = sisterList_new.distanceAligned; % copy the distances

sisterList_new = rmfield(sisterList_new,'coords1Aligned');
sisterList_new = rmfield(sisterList_new,'coords2Aligned'); % remove the alignedCoord fields
sisterList_new = rmfield(sisterList_new,'distanceAligned');


%% dependencies
% sisterList_pole_aligned depends on dataProperites, initCoords, planeFit, tracks,
% sisterList, updatedClass, frameAlignment, poles, poleReferenceFrame,
% planeFit_pole, tracks_pole, sisterList_pole, updatedClass_pole and frameAlignment_pole
% dependencies = struct('dataProperties',[],'initCoord',[],'planeFit',[],'tracks',[],'sisterList',[],'updatedClass',[],'frameAlignment',[],'poles',[],'poleReferenceFrame',[],'planeFit_pole',[],'tracks_pole',[],'sisterList_pole',[],'updatedClass_pole',[],'frameAlignment_pole',[]);

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
planeFit_poleName = dataStruct.planeFit_poleName;
planeFit_poleV = getVersion(planeFit_poleName);
tracks_poleName = dataStruct.tracks_poleName;
tracks_poleV = getVersion(tracks_poleName);
sisterList_poleName = dataStruct.sisterList_poleName;
sisterList_poleV = getVersion(sisterList_poleName);
updatedClass_poleName = dataStruct.updatedClass_poleName;
updatedClass_poleV = getVersion(updatedClass_poleName);
frameAlignment_poleName = dataStruct.frameAlignment_poleName;
frameAlignment_poleV = getVersion(frameAlignment_poleName);

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
dependencies.planeFit_pole = planeFit_poleV;
dependencies.tracks_pole = tracks_poleV;
dependencies.sisterList_pole = sisterList_poleV;
dependencies.updatedClass_pole = updatedClass_poleV;
dependencies.frameAlignment_pole = frameAlignment_poleV;

% save into the first sisterList_pole_aligned
sisterList_new(1).dependencies = dependencies;


dataStruct.sisterList_pole_aligned = sisterList_new; % save the reults

end

