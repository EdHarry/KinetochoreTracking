function dataStruct = generateSisterGroupingParam( dataStruct, useAlignment, maxAngle, maxDist, minOverlap, useAnaphase, robust)


%process user input
groupSisters.useAlignment = useAlignment;
groupSisters.maxAngle = maxAngle;
groupSisters.maxDist = maxDist;
groupSisters.minOverlap = minOverlap;
groupSisters.useAnaphase = useAnaphase;
groupSisters.robust = robust;

%save in dataStruct
dataStruct.dataProperties.groupSisters = groupSisters;

end

