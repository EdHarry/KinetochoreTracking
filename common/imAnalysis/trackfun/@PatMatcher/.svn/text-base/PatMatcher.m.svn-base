function pm=PatMatcher(tImg,sImg)
%Generate a pattern matcher object
%
%SYNOPSIS : pm=PatMatcher(tImg,sImg)
%
% INPUT:  tImg: template image structure
%         sImg: search image structure
%
% OUTPUT: PatternMatcher object 
%
% constructor:
%                  PatMatcher(tImg,sImg)
%
% public methods:
%                  pm = Match(pm,tPos,pos0,pSze,cmd)
%                  
% private methods:
%
% c: 4/1/00	dT

% check if there is data 
if(isempty(tImg) | isempty(sImg))
      pm.info.status = 'empty images';
   return;
end;

%init pm data
pm.info.status='Ok';

% default values
pm.cmd.params = [];
pm.cmd.interpol = 0;
pm.cmd.maxEx = [];
pm.cmd.prec = 0.001;
pm.cmd.maxIter=20;

% check the data type of the images
if(isa(tImg.data,'uint8'))
   pm.tempImg.data = double(tImg.data)/255;
else
   pm.tempImg=tImg;
end;

if(isa(sImg.data,'uint8'))
   pm.srchImg.data = double(sImg.data)/255;
else
   pm.srchImg=sImg;
end;

pm.template.values =[];
pm.template.mean=[];
pm.template.std=[];
pm.template.coord=[];

pm.model.coord=[];
pm.model.grad=[];
pm.model.params=[];
pm.model.pErr=1;

pm.patch=[];
pm.patCenter=0;
pm.pSze = 0;
pm.tPos = 0;
pm.pos0 = 0;

pm.interp={'*nearest','*linear','*cubic','spline'};

pm = class(pm,'PatMatcher');