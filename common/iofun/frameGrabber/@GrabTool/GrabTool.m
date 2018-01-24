function gbT=GrabTool()
%GRABTOOL create image acquisition tool object
%
% SYNOPSIS gbT=GrabTool
%
% INPUT:  None
%
% constructor:
%                  GrabTool()
%
% public methods:
%                  gbT = OpenGrabPanel(gbT)
%                  [fgTypeList,nrOfGrabbers] = InitFrameGrabbers(gbT)
%                  SaveImage(gbT,img)
% private methods:
%						 GrabPanelCB(gbT)
%                  GrabPanel(gbT)
%
% c: 24/8/99	dT

% init grabtool object

gbT.movie.data = [];
gbT.movie.stackIndx=[];
gbT.movie.stackIndxVal=0;
gbT.movie.cMap = [];

% length of movie in secs
gbT.movie.time = 10;
gbT.movie.fps = 1;
% movie repetition
gbT.movie.repeat = 1;
gbT.movie.maxImgInStack = 1000;
gbT.movie.buttonStr = {'Create Movie', 'Show Movie'};


gbT.fg.active=1;

% default shutter time in micro seconds
gbT.fg.shutterTime=5000;

% nr(x) = number of grabbers of type 'typeList(x)'
gbT.fg.nr=0;
gbT.fg.typeList=['None'];

gbT.acqDir = pwd;
gbT.acqFile = '';

gbT.data=[];
gbT.subFrame = [];
gbT.viewPanelH = [];
gbT.viewPanel2H = [];
gbT.grabPanelH=[];
gbT = class(gbT,'GrabTool');
