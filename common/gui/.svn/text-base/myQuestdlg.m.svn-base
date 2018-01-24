function ButtonName=myQuestdlg(varargin)
%QUESTDLG Question dialog box.
%  ButtonName=QUESTDLG(Question) creates a modal dialog box that 
%  automatically wraps the cell array or string (vector or matrix) 
%  Question to fit an appropriately sized window.  The name of the 
%  button that is pressed is returned in ButtonName.  The Title of 
%  the figure may be specified by adding a second string argument.  
%  Question will be interpreted as a normal string.  
%
%  SYNOPSIS ButtonName=questdlg(Question,Title,Btn1,Btn2,...,default,options)
%
%  INPUT    Question: main text in dialog box. To use TeX interpretation
%                     for the Question string, specify the option
%                     Interpreter with the value tex in the options
%           Title   : (opt) title of the dialog box {empty}
%           Btn1..n : (opt) button names {'yes','no','cancel'}
%           default : (opt) button name to be returned if a key is pressed
%           options : (opt) structure with figure property names as
%                     fieldnames and property values in the fields. Specify
%                     options.Interpreter = 'tex' to use TeX interpretation
%                     of the question string
%                     You can only change the offset of the figure, not its
%                     size with the figure property 'Position'
%
%  OUTPUT   ButtonName : name of the button that has been pressed. If the
%                        dialog box is closed without a valid selection,
%                        the return value is empty.
%
%  Example:
%
%  ButtonName=questdlg('What is your wish?', ...
%                      'Genie Question', ...
%                      'Food','Clothing','Money','Money');
%  
%  switch ButtonName,
%    case 'Food', 
%     disp('Take this BigMac Meal');
%    case 'Clothing',
%     disp('You would look so cute in this bathing suit')
%     case 'Money',
%      disp('Here''s a credit card');
%  end % switch
%
%  See also TEXTWRAP, INPUTDLG.
%
%  based on matlab's questdlg
%  modified by jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent DefaultValid;

if nargin<1
    error('Too few arguments for QUESTDLG');
else
    Question = varargin{1};
end
if ~iscell(Question),Question=cellstr(Question);end

if strcmp(Question{1},'#FigKeyPressFcn'),
    QuestFig=get(0,'CurrentFigure');
    AsciiVal= abs(get(QuestFig,'CurrentCharacter'));
    if ~isempty(AsciiVal),
        if AsciiVal==32 | AsciiVal==13,
            % Check if the default string matches any button string.
            % If not then dont resune till the user selects a valid input.
            if(~DefaultValid)
	        warnstate = warning('backtrace','off');
                warning('MATLAB:QUESTDLG:stringMismatch','Default string does not match any button string name.');
		warning(warnstate);
                return;
            end
            set(QuestFig,'UserData',1);
            uiresume(QuestFig);
        end %if AsciiVal
    end %if ~isempty
    return
end
%%%%%%%%%%%%%%%%%%%%%
%%% General Info. %%%
%%%%%%%%%%%%%%%%%%%%%
Black      =[0       0        0      ]/255;
LightGray  =[192     192      192    ]/255;
LightGray2 =[160     160      164    ]/255;
MediumGray =[128     128      128    ]/255;
White      =[255     255      255    ]/255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%% Nargin Check & assign defaults %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set defaults
Interpreter='none';
Title = '';
ButtonList{1} = 'yes';
ButtonList{2} = 'no';
ButtonList{3} = 'cancel';
ButtonList{4} = 'yes'; %will become Default
Default = '';

%check nargin
if nargout>1,error('Wrong number of output arguments for QUESTDLG');end

%find title
if nargin > 1 & ~isempty(varargin{2}) & ~isstruct(varargin{2})
    Title = varargin{2};
end

%find options
if isstruct(varargin{end})
    options = varargin{end};
    nbMax = nargin - 3;
    %test for field interpreter
    if isfield(options,'Interpreter')
        Interpreter = options.Interpreter;
        rmfield(options,'Interpreter'); %we can't use this field
    end
else
    options = [];
    nbMax = nargin - 2;
end


%look for buttons
if nbMax > 0
    ButtonList = varargin(3:2+nbMax);
end

%look for default
if strmatch(ButtonList{end},ButtonList(1:end-1))
    Default = ButtonList{end};
    ButtonList = ButtonList(1:end-1);
end

NumButtons = length(ButtonList);


%%%%%%%%%%%%%%%%%%%%%%%
%%% Create QuestFig %%%
%%%%%%%%%%%%%%%%%%%%%%%
FigPos=get(0,'DefaultFigurePosition');
FigWidth=190;FigHeight=50;
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','points');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);

    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','points');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
QuestFig=dialog(                                               ...
               'Visible'         ,'off'                      , ...
               'Name'            ,Title                      , ...
               'Pointer'         ,'arrow'                    , ...
               'Units'           ,'points'                   , ...
               'Position'        ,FigPos                     , ...
               'KeyPressFcn'     ,'myQuestdlg #FigKeyPressFcn;', ...
               'UserData'        ,0                          , ...
               'IntegerHandle'   ,'off'                      , ...
               'WindowStyle'     ,'normal'                   , ... 
               'HandleVisibility','callback'                 , ...
               'Tag'             ,Title                        ...
               );

%%%%%%%%%%%%%%%%%%%%%
%%% Set Positions %%%
%%%%%%%%%%%%%%%%%%%%%
%DefOffset=3;
DefOffset=7;

%IconWidth=32;
IconWidth=36;
%IconHeight=32;
IconHeight=37;
IconXOffset=DefOffset;
IconYOffset=FigHeight-DefOffset-IconHeight;
IconCMap=[Black;get(QuestFig,'Color')];

DefBtnWidth=40;
BtnHeight=20;
BtnYOffset=DefOffset;
BtnFontSize=get(0,'FactoryUIControlFontSize');
BtnFontName=get(0,'FactoryUIControlFontName');

BtnWidth=DefBtnWidth;

ExtControl=uicontrol(QuestFig   , ...
                     'Style'    ,'pushbutton', ...
                     'String'   ,' '         , ...
                     'FontUnits','points'    , ...                     
                     'FontSize' ,BtnFontSize , ...
                     'FontName' ,BtnFontName   ...
                     );

%set button width                     
for nBtn = 1:NumButtons                     
    set(ExtControl,'String',ButtonList{nBtn});
    BtnExtent=get(ExtControl,'Extent');
    BtnWidth=max(BtnWidth,BtnExtent(3)+8);
end

delete(ExtControl);

MsgTxtXOffset=IconXOffset+IconWidth;

FigWidth=max(FigWidth,MsgTxtXOffset+NumButtons*(BtnWidth+2*DefOffset));
FigPos(3)=FigWidth;
set(QuestFig,'Position',FigPos);

BtnXOffset=zeros(NumButtons,1);

% distribute the buttons
% IXOff-AddOff-B1-DefOff-AddOff-B2-...B3-AddOff-IXOff

%calculate AddOffset
AddOffset = max((FigWidth-(BtnWidth+DefOffset)*NumButtons+DefOffset-2*IconXOffset)/(NumButtons+1),0);

%distribute: IXOffset-Btn1-DefOffset-AddOffset-Btn2-...-BtnX-IXOffset
for nBtn = 1:NumButtons
    %Offset          = Left Margin           + Width of all buttons & margins to the left
    BtnXOffset(nBtn) = IconXOffset + AddOffset +(nBtn-1)*(BtnWidth+DefOffset+AddOffset);
end

MsgTxtYOffset=DefOffset+BtnYOffset+BtnHeight;
MsgTxtWidth=FigWidth-DefOffset-MsgTxtXOffset-IconWidth;
MsgTxtHeight=FigHeight-DefOffset-MsgTxtYOffset;
MsgTxtForeClr=Black;
MsgTxtBackClr=get(QuestFig,'Color');

CBString='uiresume(gcf)';

% Checks to see if the Default string passed does match one of the
% strings on the buttons in the dialog. 
DefaultValid = 0;
for nBtn = 1:NumButtons
    ButtonString=ButtonList{nBtn};
    ButtonTag=['Btn',num2str(nBtn)];
    BtnHandle(nBtn)=uicontrol(QuestFig            , ...
			'Style'              ,'pushbutton', ...
			'Units'              ,'points'    , ...
			'Position'           ,[ BtnXOffset(nBtn) BtnYOffset  ...
		                                BtnWidth       BtnHeight   ...
		                              ]           , ...
			'CallBack'           ,CBString    , ...
			'String'             ,ButtonString, ...
			'HorizontalAlignment','center'    , ...
			'FontUnits'          ,'points'    , ...
			'FontSize'           ,BtnFontSize , ...
			'FontName'           ,BtnFontName , ...
			'Tag'                ,ButtonTag     ...
			);
    if ~isempty(Default) & strcmp(ButtonString, Default)
        DefaultValid = 1;
        set(BtnHandle(nBtn),'FontWeight','bold');
    end
end



MsgHandle=uicontrol(QuestFig            , ...
                   'Style'              ,'text'         , ...
                   'Units'              ,'points'       , ...
                   'Position'           ,[MsgTxtXOffset      ...
                                          MsgTxtYOffset      ...
                                          0.95*MsgTxtWidth   ...
                                          MsgTxtHeight       ...
                                         ]              , ...
                   'String'             ,{' '}          , ...
                   'Tag'                ,'Question'     , ...
                   'HorizontalAlignment','left'         , ...    
                   'FontUnits'          ,'points'       , ...
                   'FontWeight'         ,'bold'         , ...
                   'FontSize'           ,BtnFontSize    , ...
                   'FontName'           ,BtnFontName    , ...
                   'BackgroundColor'    ,MsgTxtBackClr  , ...
                   'ForegroundColor'    ,MsgTxtForeClr    ...
                   );

[WrapString,NewMsgTxtPos]=textwrap(MsgHandle,Question,75);

NumLines=size(WrapString,1);

% The +2 is to add some slop for the border of the control.
MsgTxtWidth=max(MsgTxtWidth,NewMsgTxtPos(3)+2);
MsgTxtHeight=NewMsgTxtPos(4)+2;

MsgTxtXOffset=IconXOffset+IconWidth+DefOffset;
%adjust figureWidth
FigWidth=max((BtnWidth+DefOffset)*NumButtons-DefOffset+2*IconXOffset, ...
             MsgTxtXOffset+MsgTxtWidth+DefOffset);

        
% Center Vertically around icon  
if IconHeight>MsgTxtHeight,
  IconYOffset=BtnYOffset+BtnHeight+DefOffset;
  MsgTxtYOffset=IconYOffset+(IconHeight-MsgTxtHeight)/2;
  FigHeight=IconYOffset+IconHeight+DefOffset;    
% center around text    
else,
  MsgTxtYOffset=BtnYOffset+BtnHeight+DefOffset;
  IconYOffset=MsgTxtYOffset+(MsgTxtHeight-IconHeight)/2;
  FigHeight=MsgTxtYOffset+MsgTxtHeight+DefOffset;    
end    

%calculate AddOffset again, because the real text has now been added
AddOffset = max((FigWidth-(BtnWidth+DefOffset)*NumButtons+DefOffset-2*IconXOffset)/(NumButtons+1),0);

%distribute: IXOffset-Btn1-DefOffset-AddOffset-Btn2-...-BtnX-IXOffset
for nBtn = 1:NumButtons
    %Offset          = Left Margin           + Width of all buttons & margins to the left
    BtnXOffset(nBtn) = IconXOffset + AddOffset +(nBtn-1)*(BtnWidth+DefOffset+AddOffset);
end


FigPos(3:4)=[FigWidth FigHeight];

set(QuestFig ,'Position',FigPos);

BtnPos=get(BtnHandle,{'Position'});BtnPos=cat(1,BtnPos{:});
BtnPos(:,1)=BtnXOffset;
BtnPos=num2cell(BtnPos,2);  
set(BtnHandle,{'Position'},BtnPos);  

delete(MsgHandle);
AxesHandle=axes('Parent',QuestFig,'Position',[0 0 1 1],'Visible','off');

MsgHandle=text( ...
    'Parent'              ,AxesHandle                      , ...
    'Units'               ,'points'                        , ...
    'FontUnits'           ,'points'                        , ...
    'FontSize'            ,BtnFontSize                     , ...
    'FontName'            ,BtnFontName                     , ...
    'HorizontalAlignment' ,'left'                          , ...
    'VerticalAlignment'   ,'bottom'                        , ...
    'HandleVisibility'    ,'callback'                      , ...
    'Position'            ,[MsgTxtXOffset MsgTxtYOffset 0] , ...
    'String'              ,WrapString                      , ...
    'Interpreter'         ,Interpreter                     , ...
    'Tag'                 ,'Question'                        ...
    );

IconAxes=axes(                                      ...
             'Units'       ,'points'              , ...
             'Parent'      ,QuestFig              , ...  
             'Position'    ,[IconXOffset IconYOffset  ...
                             IconWidth IconHeight], ...
             'NextPlot'    ,'replace'             , ...
             'Tag'         ,'IconAxes'              ...
             );         
 
set(QuestFig ,'NextPlot','add');

%set options
if ~isempty(options)
    propertyNames = fieldnames(options);
    propertyValues = struct2cell(options);
    for nopt = 1:length(options)
        try
            if strcmp(propertyNames{nopt},'Position')
                qPos = get(QuestFig,'Position');
                newPos = propertyValues{nopt};
                qPos(1:2) = newPos(1:2);
                set(QuestFig,'Position',qPos);
            end
            set(QuestFig,propertyNames{nopt},propertyValues{nopt});
        catch
            disp('Warning: ',propertyNames{nopt},' did is either not a valid figure option or contains an invalid value!')
        end
    end
end

load dialogicons.mat
IconData=questIconData;
questIconMap(256,:)=get(QuestFig,'color');
IconCMap=questIconMap;

Img=image('CData',IconData,'Parent',IconAxes);
set(QuestFig, 'Colormap', IconCMap);
set(IconAxes, ...
   'Visible','off'           , ...
   'YDir'   ,'reverse'       , ...
   'XLim'   ,get(Img,'XData'), ...
   'YLim'   ,get(Img,'YData')  ...
   );
set(findobj(QuestFig),'HandleVisibility','callback');
set(QuestFig ,'WindowStyle','modal','Visible','on');



drawnow;

%wait for user input
uiwait(QuestFig);

TempHide=get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');

if any(get(0,'Children')==QuestFig),
  if get(QuestFig,'UserData'),
    ButtonName=Default;
  else,
    ButtonName=get(get(QuestFig,'CurrentObject'),'String');
  end
  delete(QuestFig);
else
  ButtonName='';
end

set(0,'ShowHiddenHandles',TempHide);
