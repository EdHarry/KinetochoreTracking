function uiProfDeleteFcn
% ui delete lines from original picture

proFig=gcbf;
lineH=get(proFig,'UserData');
delete(lineH);