function aImarisApplication = StartImaris_v73
% EHarry Nov 2011
javaaddpath '/Applications/Imaris 7.3.0.app/Contents/SharedSupport/XT/matlab/ImarisLib.jar';
vImarisLib = ImarisLib;
vObjectId = StartImarisInstance_v73;
for vIndex = 1:1000 % do several attempts waiting for Imaris to be registered
    aImarisApplication = vImarisLib.GetApplication(vObjectId);
    if ~isempty(aImarisApplication)
        break
    end
    pause(0.1)
end
end