function modifyTiffOptions(adepth,res,msg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  SETTING VIEW PANEL TIFF OPTIONS AND SAVE MENU 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Activates 'x bit' menu in VIEWPANEL and 'save stack' menu in TRACK PANEL if normalization successful

% All 'x bit' options off
if res~=2
   par=findobj('Tag','UIVIEWMENU_FILE_OPTIONS');
   ch=get(par,'Children');
   
   for counter1=1:length(ch)
      set(ch(counter1),'Checked','Off');
   end 
end

% What to turn on depends on 'res'
if res==0
   % Activates current 'x bit' menu option
       
   switch adepth
   case 8
      counter1=9;
   case 9
      counter1=8;
   case 10
      counter1=7;
   case 11
      counter1=6;
   case 12
      counter1=5;
   case 13
      counter1=4;
   case 14
      counter1=3;
   case 15
      counter1=2;
   case 16
      counter1=1;
   otherwise
      % If the image is not in the range 8..16 bit the menu is set to 'Auto search (smart)'
      counter1=10;
   end
   set(ch(counter1),'Checked','On');                                                                                                                                                                                                         
   
   % Activates Save Stack menu  if normalization successful
   par=findobj('Tag','UIVIEWMENU_FILE_SAVE');
   ch=get(par,'Children');
   for counter1=1:length(ch)
      set(ch(counter1),'Enable','On');
   end  
end

if res==-1
   % 'x bit' set to 'auto search (smart)'
   counter1=10;
   set(ch(counter1),'Checked','On');
   
   % save disabled
   par=findobj('Tag','UISSDTRACKMENU_FILE_SAVESTACK');
   set(par,'Enable','Off');
end                                                                                                                                                                                                                                       

% Turns on 'Save' and 'Save Stack'
set(findobj(gcbf,'Tag','UISSDTRACKMENU_FILE_SAVE'),'Enable','On');
set(findobj(gcbf,'Tag','UISSDTRACKMENU_FILE_SAVESTACK'),'Enable','On');

% Returns msg from imreadnd
if ~isempty(msg)
   displ(msg);
end
