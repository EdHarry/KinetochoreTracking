function aServer = GetServer
% EHarry Nov 2011
javaaddpath '/Applications/Imaris 7.3.0.app/Contents/SharedSupport/XT/matlab/ImarisLib.jar';
vImarisLib = ImarisLib;
aServer = vImarisLib.GetServer;
end