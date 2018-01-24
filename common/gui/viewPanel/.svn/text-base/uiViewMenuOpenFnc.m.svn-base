function uiViewMenuOpenFnc
% callback function
%
% m: 17/6/99	dT

% Reads Tiff Options
[adepth,algorithm,res,msg]=readTiffOptions;

% Loading and normalizing image
[fName,dirName] = uigetfile('*.tif','View Panel: Open ...');
if(isa(fName,'char') & isa(dirName,'char'))
   cd(dirName);
   % Display
   displ('Evaluating... please wait');
   % Starting normalization
   [I,adepth,msg,res]=imreadnd([dirName,fName],adepth,algorithm);
   uiViewPanelShowImg(I,0,gcbf);
end;

% Modifies Tiff options and save menu
modifyTiffOptions(adepth,res,msg);
