function uiViewMenuScaleFnc
%
% callback for scale menu
figH = gcbf;

% check whether there exists a scale
myScale = GetUserData(figH,'myScale');
if (isempty(myScale))
   myScale=Scale(figH);
end;
SetScale(myScale);
