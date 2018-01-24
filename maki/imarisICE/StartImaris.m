function aImarisApplication = StartImaris
% EHarry Nov 2011
javaaddpath '/Applications/Imaris 7.4.0.app/Contents/SharedSupport/XT/matlab/ImarisLib.jar';
vImarisLib = ImarisLib;
vObjectId = StartImarisInstance;
for vIndex = 1:1000 % do several attempts waiting for Imaris to be registered
    aImarisApplication = vImarisLib.GetApplication(vObjectId);
    if ~isempty(aImarisApplication)
        break
    end
    pause(0.1)
end
end