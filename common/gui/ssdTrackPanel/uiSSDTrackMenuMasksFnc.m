function uiSSDTrackMenuMasksFnc
%UISSDTRACKMENUMASKSFNC callback for the mask menu of the SSD track panel

whoCalls = get(gcbo,'Tag');

switch whoCalls
case 'UISSDTRACKMENU_MASK', masterMenu;
case 'UISSDTRACKMENU_MASK_NO', maskChange(0);
case 'UISSDTRACKMENU_MASK_MANCTR', maskChange(1);
case 'UISSDTRACKMENU_MASK_MANSEGM', maskChange(2);
case 'UISSDTRACKMENU_MASK_DICLF', maskChange(3);
case 'UISSDTRACKMENU_MASK_DICAD', maskChange(4);
case 'UISSDTRACKMENU_MASK_POSLF', maskChange(5);
case 'UISSDTRACKMENU_MASK_POSAD', maskChange(6);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions

function masterMenu

pI = get(gcbf,'UserData');

% uncheck the currently checked menu point
set(findobj(gcbo,'Checked','on'),'Checked','off');

switch pI.tImgMaskType
case 0, set(findobj(gcbo,'Tag','UISSDTRACKMENU_MASK_NO'),'Checked','on');
case 1, set(findobj(gcbo,'Tag','UISSDTRACKMENU_MASK_MANCTR'),'Checked','on');
case 2, set(findobj(gcbo,'Tag','UISSDTRACKMENU_MASK_MANSEGM'),'Checked','on');
case 3, set(findobj(gcbo,'Tag','UISSDTRACKMENU_MASK_DICLF'),'Checked','on');
case 4, set(findobj(gcbo,'Tag','UISSDTRACKMENU_MASK_DICAD'),'Checked','on');
case 5, set(findobj(gcbo,'Tag','UISSDTRACKMENU_MASK_POSLF'),'Checked','on');
case 6, set(findobj(gcbo,'Tag','UISSDTRACKMENU_MASK_POSAD'),'Checked','on');
otherwise, error('invalid mask type stored');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maskChange(code)

pI = get(gcbf,'UserData');

pI.tImgMaskType = code;
set(gcbf,'UserData',pI);

if(~isempty(pI.tImg))
   msgbox('Changes in mask setting will be activated after the next reinitialization of template',...
      'SSD Track Message...','modal');
end;
