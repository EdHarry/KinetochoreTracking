function aImarisApplication = GetImaris
% EHarry Nov 2011
javaaddpath '/Applications/Imaris 7.3.0.app/Contents/SharedSupport/XT/matlab/ImarisLib.jar';
vImarisLib = ImarisLib;
vObjectId = 0; % this might be replaced by "vObjectId = <a name=getobjectid><b>GetObjectId</b></a>" (see later)
aImarisApplication = vImarisLib.GetApplication(vObjectId);
end