function dataStruct = makiUpdateClassPole_oldVersion( dataStruct )
%MAKIUPDATECLASSPOLE updates the class and phase of a cell ignoring the
%poles, which otherwise would incorrectly be identified as (outling) kinetochores
% EHarry Jan 2012

dataStruct_temp = dataStruct; % copy the dataStruct for the operations

dataStruct_temp.planeFit = []; % delete the original planeFit

dataStruct_temp = makiConvertPoleCoordsToInitCoords(dataStruct_temp); % convert the poleCoors to initCoords

try
    dataStruct_temp = makiFitPlane_CENPQ(dataStruct_temp,0); % try to fit a new plane
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


dataStruct.planeFit_pole = planeFit; % save as a new struct

end

