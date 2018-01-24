function [xi,yi] = imSelCurve(varargin)
%imSelCurve: Interactively select a curve in the current figure.
%
%    Select a curve on the current axes by clicking points in the figure. 
%    Exit when you hit the return key.
%
% SYNTAX: 
%    [xi,yi] = imSelCurve
%       The coordinates of the curve is returned in 'xi' and 'yi'.
%    [xi,yi] = imSelCurve(maxNumPoints)
%       Specify the maximum number of break points on the curve. When this
%       number is reached by mouse clicking, the program exits and return the
%       curve. 'maxNumPoints' has to be 0 or positive. When it is zero, you
%       can click as many points as possible and the program exits when you
%       hit the return key. 

if nargin > 1
   error('Too many input arguments.');
end

if nargin == 1
   if varargin{1} < 0
      error(['The maximum number of specified break points has to ' ...
         'be zero or positive.']);
   end
   handles.maxNumPoints = varargin{1};
else
   handles.maxNumPoints = 0;
end

%Get the current figure.
h = gcf; figure(h);

%Save the old guidata and some properties of the current figure.
oldHandles         = guidata(h);
oldPointer         = get(h,'Pointer');
oldButtonDownFcn   = get(h,'WindowButtonDownFcn');
oldButtonMotionFcn = get(h,'WindowButtonMotionFcn');
oldKeyPressFcn     = get(h,'KeyPressFcn');

%Initialize with empty data.
handles.numPoints = 0;
handles.xi        = [];
handles.yi        = [];
handles.dashL     = {};
handles.dotL      = {};
handles.mDashL    = [];
handles.mDotL     = [];

%Update GUI data.
guidata(h,handles);

set(h,'DoubleBuffer','on','Pointer','cross');
set(h,'WindowButtonDownFcn',@buttonDown_Callback);
set(h,'WindowButtonMotionFcn',@buttonMove_Callback);
set(h,'KeyPressFcn',@keyPress_Callback);

%Set it to be modal.
%set(h,'WindowStyle','modal','MenuBar','figure');

w       = waitforbuttonpress;
ckey    = get(h,'CurrentCharacter');
handles = guidata(h);
while (w == 0 & (handles.maxNumPoints == 0 | handles.numPoints < ...
   handles.maxNumPoints)) | (w == 1 & double(ckey) ~= 13)
   %When there is no key stroke from the key board or the key stroke is not
   % return, keep waiting.
   w       = waitforbuttonpress;
   ckey    = get(h,'CurrentCharacter');
   handles = guidata(h);
end

%Block mouse action and key pressing call back.
set(h,'WindowStyle','normal','Pointer',oldPointer);
set(h,'WindowButtonDownFcn','');
set(h,'WindowButtonMotionFcn','');
set(h,'KeyPressFcn','');

handles = guidata(h);

xi = handles.xi;
yi = handles.yi;

%Delete the lines that are drawn during selection.
if handles.numPoints > 1
   for k = 1:handles.numPoints-1
      if ishandle(handles.dashL{k})
         delete(handles.dashL{k});
      end
      if ishandle(handles.dotL{k})
         delete(handles.dotL{k});
      end
   end
end

if ishandle(handles.mDashL)
   delete(handles.mDashL);
   handles.mDashL = [];
end
if ishandle(handles.mDotL)
   delete(handles.mDotL);
   handles.mDotL = [];
end

%Set guidata and figure properties back to the state before calling this
% function.
handles = oldHandles;
guidata(h,handles);
set(h,'WindowStyle','normal','Pointer',oldPointer);
set(h,'WindowButtonDownFcn',oldButtonDownFcn);
set(h,'WindowButtonMotionFcn',oldButtonMotionFcn);
set(h,'KeyPressFcn',oldKeyPressFcn);


% --- Call back function for WindowButtonDownFcn.
function buttonDown_Callback(obj,eventdata)

%Get the GUI data
handles = guidata(obj);

numPoints = handles.numPoints+1;

%Delete previously drawn line while mounse is moving.
if ishandle(handles.mDashL)
   delete(handles.mDashL);
end
if ishandle(handles.mDotL)
   delete(handles.mDotL);
end
handles.mDashL = [];
handles.mDotL  = [];

p = get(gca,'CurrentPoint');

if numPoints > 1
   handles.dashL{numPoints-1} = line([handles.xi(end);p(1,1)], ...
      [handles.yi(end);p(1,2)]);
   set(handles.dashL{numPoints-1},'LineStyle','--');
   set(handles.dashL{numPoints-1},'color','k');
   handles.dotL{numPoints-1} = line([handles.xi(end);p(1,1)], ...
      [handles.yi(end);p(1,2)]);
   set(handles.dotL{numPoints-1},'LineStyle',':');
   set(handles.dotL{numPoints-1},'color','w');
end

handles.xi(numPoints) = p(1,1);
handles.yi(numPoints) = p(1,2);

handles.numPoints = numPoints;

%Update GUI data.
guidata(obj,handles);

% --- Call back function for WindowButtonMotionFcn.
function buttonMove_Callback(obj,eventdata)

%Get the GUI data
handles = guidata(obj);

if handles.numPoints == 0
   return;
end

%Delete previously drawn line while mounse is moving.
if ishandle(handles.mDashL)
   delete(handles.mDashL);
end
if ishandle(handles.mDotL)
   delete(handles.mDotL);
end
handles.mDashL = [];
handles.mDotL  = [];

p = get(gca,'CurrentPoint');

%Draw the moving line
handles.mDashL = line([handles.xi(end);p(1,1)],[handles.yi(end);p(1,2)]);
set(handles.mDashL,'LineStyle','--');
set(handles.mDashL,'color','k');
handles.mDotL = line([handles.xi(end);p(1,1)],[handles.yi(end);p(1,2)]);
set(handles.mDotL,'LineStyle',':');
set(handles.mDotL,'color','w');

%Update GUI data.
guidata(obj,handles);

% --- Call back function for KeyPressFcn.
function keyPress_Callback(obj,eventdata)

handles = guidata(obj);

numPoints = handles.numPoints;
if numPoints == 0
   return;
end

p    = get(gca,'CurrentPoint');
ckey = get(gcf,'CurrentCharacter');

if double(ckey) == 8
   %When the key stroke is backspace (ascii code 8), delete the last selected 
   % point.
   handles.xi(end)   = [];
   handles.yi(end)   = [];
   handles.numPoints = numPoints-1;

   %Delete the last line.
   if numPoints > 1
      if ishandle(handles.dashL{end})
         delete(handles.dashL{end});
      end
      if ishandle(handles.dotL{end})
         delete(handles.dotL{end});
      end
      handles.dashL(end) = [];
      handles.dotL(end)  = [];
   end

   if ishandle(handles.mDashL)
      delete(handles.mDashL);
      handles.mDashL = [];
   end
   if ishandle(handles.mDotL)
      delete(handles.mDotL);
      handles.mDotL  = [];
   end

   if numPoints > 1
      handles.mDashL = line([handles.xi(end);p(1,1)],[handles.yi(end);p(1,2)]);
      set(handles.mDashL,'LineStyle','--');
      set(handles.mDashL,'color','k');
      handles.mDotL = line([handles.xi(end);p(1,1)],[handles.yi(end);p(1,2)]);
      set(handles.mDotL,'LineStyle',':');
      set(handles.mDotL,'color','w');
   end
end

%Update GUI data.
guidata(obj,handles);

