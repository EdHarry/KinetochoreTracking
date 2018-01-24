function aObjectId = GetObjectId
% EHarry Nov 2011
vServer = GetServer;
vNumberOfObjects = vServer.GetNumberOfObjects;
for vIndex = 0:vNumberOfObjects-1
    aObjectId = vServer.GetObjectID(vIndex);
    break % work with the ID, return first one
end
end