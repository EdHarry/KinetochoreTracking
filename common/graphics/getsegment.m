function [x,y] = getsegment(varargin)
%GETSEGMENT Select a straight line segment with mouse.
%
% SYNOPSIS [x,y] = getsegment(varargin)
%
% INPUT
%   varargin:  [X,Y] = getsegment(FIG) lets you select a segment in the
%                      current axes of figure FIG using the mouse.   
%                      Use normal button clicks to add points to the segement.
%                      Pressing BACKSPACE or DELETE removes the previously 
%                      selected point from the polyline.
%
%              [X,Y] = getsegment(AX) lets you select a segment in the axes
%                       specified by the handle AX.
%
%              [X,Y] = getsegment is the same as [X,Y] = getsegment(GCF).
%
%  OUTPUT
%    Coordinates of the segement are returned in X and Y. 
%
%  NOTE function is a modification of getpts() in the image toolbox

% modified getline code, which is a native MATLAB image 
% toolbox function
%   Callback syntaxes:
%        getsegement('KeyPress')
%        getsegment('FirstButtonDown')
%        getsegment('SecondButtonDown')
%        getsegment('ButtonMotion')


global GETLINE_FIG GETLINE_AX GETLINE_H1 GETLINE_H2
global GETLINE_X GETLINE_Y


if ((length(varargin) >= 1) & isstr(varargin{1}))
    % Callback invocation: 'KeyPress', 'FirstButtonDown',
    % 'SecondButtonDown', or 'ButtonMotion'.
    feval(varargin{:});
    return;
end

GETLINE_X = [];
GETLINE_Y = [];

if (length(varargin) < 1)
    GETLINE_AX = gca;
    GETLINE_FIG = get(GETLINE_AX, 'Parent');
else
    if (~ishandle(varargin{1}))
        error('First argument is not a valid handle');
    end
    
    switch get(varargin{1}, 'Type')
    case 'figure'
        GETLINE_FIG = varargin{1};
        GETLINE_AX = get(GETLINE_FIG, 'CurrentAxes');
        if (isempty(GETLINE_AX))
            GETLINE_AX = axes('Parent', GETLINE_FIG);
        end

    case 'axes'
        GETLINE_AX = varargin{1};
        GETLINE_FIG = get(GETLINE_AX, 'Parent');

    otherwise
        error('First argument should be a figure or axes handle');

    end
end

% Bring target figure forward
figure(GETLINE_FIG);

% Remember initial figure state
state= uisuspend(GETLINE_FIG);

% Set up initial callbacks for initial stage
set(GETLINE_FIG, 'Pointer', 'crosshair');
set(GETLINE_FIG, 'WindowButtonDownFcn', 'getsegment(''FirstButtonDown'');');
set(GETLINE_FIG, 'KeyPressFcn', 'getsegment(''KeyPress'');');

% Initialize the lines to be used for the drag
GETLINE_H1 = line('Parent', GETLINE_AX, ...
                  'XData', GETLINE_X, ...
                  'YData', GETLINE_Y, ...
                  'Visible', 'off', ...
                  'Clipping', 'off', ...
                  'Color', 'k', ...
                  'LineStyle', '-', ...
                  'EraseMode', 'xor');

GETLINE_H2 = line('Parent', GETLINE_AX, ...
                  'XData', GETLINE_X, ...
                  'YData', GETLINE_Y, ...
                  'Visible', 'off', ...
                  'Clipping', 'off', ...
                  'Color', 'w', ...
                  'LineStyle', '--', ...
                  'EraseMode', 'xor');

% We're ready; wait for the user to do the drag
% Wrap the call to waitfor in try-catch so we'll
% have a chance to clean up after ourselves.
errCatch = 0;
try
    waitfor(GETLINE_H1, 'UserData', 'Completed');
catch
    errCatch = 1;
end

% After the waitfor, if GETLINE_H1 is still valid
% and its UserData is 'Completed', then the user
% completed the drag.  If not, the user interrupted
% the action somehow, perhaps by a Ctrl-C in the
% command window or by closing the figure.

if (errCatch == 1)
    errStatus = 'trap';
    
elseif (~ishandle(GETLINE_H1) | ...
            ~strcmp(get(GETLINE_H1, 'UserData'), 'Completed'))
    errStatus = 'unknown';
    
else
    errStatus = 'ok';
    x = GETLINE_X(:);
    y = GETLINE_Y(:);
    % If no points were selected, return rectangular empties.
    % This makes it easier to handle degenerate cases in
    % functions that call getline.
    if (isempty(x))
        x = zeros(0,1);
    end
    if (isempty(y))
        y = zeros(0,1);
    end
end

% Delete the animation objects
if (ishandle(GETLINE_H1))
    delete(GETLINE_H1);
end
if (ishandle(GETLINE_H2))
    delete(GETLINE_H2);
end

% Restore the figure's initial state
if (ishandle(GETLINE_FIG))
   uirestore(state);
end

% Clean up the global workspace
clear global GETLINE_FIG GETLINE_AX GETLINE_H1 GETLINE_H2
clear global GETLINE_X GETLINE_Y

% Depending on the error status, return the answer or generate
% an error message.
switch errStatus
case 'ok'
    % nothing to do    
case 'trap'
    % An error was trapped during the waitfor
    error('Interruption during mouse selection.');
    
case 'unknown'
    % User did something to cause the polyline selection to
    % terminate abnormally.  For example, we would get here
    % if the user closed the figure in the middle of the selection.
    error('Interruption during mouse selection.');
end

%--------------------------------------------------
% Subfunction KeyPress
%--------------------------------------------------
function KeyPress

global GETLINE_FIG GETLINE_AX GETLINE_H1 GETLINE_H2
global GETLINE_PT1 
global GETLINE_X GETLINE_Y

key = real(get(GETLINE_FIG, 'CurrentCharacter'));
switch key
case {8, 127}  % delete and backspace keys
    % remove the previously selected point
    switch length(GETLINE_X)
    case 0
        % nothing to do
    case 1
        GETLINE_X = [];
        GETLINE_Y = [];
        % remove point and start over
        set([GETLINE_H1 GETLINE_H2], ...
                'XData', GETLINE_X, ...
                'YData', GETLINE_Y);
        set(GETLINE_FIG, 'WindowButtonDownFcn', ...
                'getsegment(''FirstButtonDown'');', ...
                'WindowButtonMotionFcn', '');
    otherwise
        % remove last point
        GETLINE_X(end) = [];
        GETLINE_Y(end) = [];
        
        set([GETLINE_H1 GETLINE_H2], ...
                'XData', GETLINE_X, ...
                'YData', GETLINE_Y);
    end
end

%--------------------------------------------------
% Subfunction FirstButtonDown
%--------------------------------------------------
function FirstButtonDown

global GETLINE_FIG GETLINE_AX GETLINE_H1 GETLINE_H2
global GETLINE_X GETLINE_Y


[GETLINE_X, GETLINE_Y] = getcurpt(GETLINE_AX);
set([GETLINE_H1 GETLINE_H2], ...
        'XData', GETLINE_X, ...
        'YData', GETLINE_Y, ...
        'Visible', 'on');

if (~strcmp(get(GETLINE_FIG, 'SelectionType'), 'normal'))
    % We're done!
    set(GETLINE_H1, 'UserData', 'Completed');
else
    % Let the motion functions take over.
    set(GETLINE_FIG, 'WindowButtonMotionFcn', 'getsegment(''ButtonMotion'');', ...
            'WindowButtonDownFcn', 'getsegment(''SecondButtonDown'');');
end

%--------------------------------------------------
% Subfunction SecondButtonDown
%--------------------------------------------------
function SecondButtonDown

global GETLINE_FIG GETLINE_AX GETLINE_H1 GETLINE_H2
global GETLINE_X GETLINE_Y

selectionType = get(GETLINE_FIG, 'SelectionType');
if (~strcmp(selectionType, 'open'))
    % We don't want to add a point on the second click
    % of a double-click

    [x,y] = getcurpt(GETLINE_AX);
    GETLINE_X = [GETLINE_X x];
    GETLINE_Y = [GETLINE_Y y];
    
    set([GETLINE_H1 GETLINE_H2], 'XData', GETLINE_X, ...
       'YData', GETLINE_Y);
    
    set(GETLINE_H1, 'UserData', 'Completed');

end

if (~strcmp(get(GETLINE_FIG, 'SelectionType'), 'normal'))
    % We're done!
    set(GETLINE_H1, 'UserData', 'Completed');
end

%-------------------------------------------------
% Subfunction ButtonMotion
%-------------------------------------------------
function ButtonMotion

global GETLINE_FIG GETLINE_AX GETLINE_H1 GETLINE_H2
global GETLINE_X GETLINE_Y

[newx, newy] = getcurpt(GETLINE_AX);

x = [GETLINE_X newx];
y = [GETLINE_Y newy];

set([GETLINE_H1 GETLINE_H2], 'XData', x, 'YData', y);

