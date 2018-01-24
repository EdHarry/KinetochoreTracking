function dataStruct = makiMakeAlignedSisterList_plate2pole( dataStruct )
%MAKIMAKEALIGNEDSISTERLIST Makes a sisterList_pole based only on alignedCoords from a pole based system and adds it to the dataStruct
% EHarry Jan 2012

try
   sisterList_new = makiConstructAlignedSisters_plate2pole(dataStruct); % make the aligned sisterList 
catch
    dataStruct.sisterList_plate2pole = [];
    return
end

if isempty(sisterList_new)
    dataStruct.sisterList_plate2pole = [];
    return
end

[sisterList_new.coords1] = sisterList_new.coords1Aligned; % copy the aligned coords to the regualar coords
[sisterList_new.coords2] = sisterList_new.coords2Aligned;

[sisterList_new.distances] = sisterList_new.distanceAligned; % copy the distances

sisterList_new = rmfield(sisterList_new,'coords1Aligned');
sisterList_new = rmfield(sisterList_new,'coords2Aligned'); % remove the alignedCoord fields
sisterList_new = rmfield(sisterList_new,'distanceAligned');

dataStruct.sisterList_plate2pole = sisterList_new; % save the reults

end

