function pt = getpt(varargin)
%GETPT Select a single point with mouse.
%
% SYNOPSIS pt = getpt(varargin)
% 
% INPUT varargin : getpt(FIG) lets you choose a points in the
%                  current axes of figure FIG using the mouse. 
%
%                  getpt(AX) lets you choose a point in the axes
%                  specified by the handle AX.
%
%                  getpt is the same as getpt(GCF).
%
% OUTPUT pt : coordinates of the selcted point
%
% NOTE this function works only if an image is displaid in the figure
%      (function is a modification of getpts() in the image toolbox)

% modified getpts code, which is a native MATLAB image 
% toolbox function
%   Callback syntaxes:
%       getpt('ButtonDown')


global GETPT_FIG GETPT_AX GETPT_H1
global GETPT_PT1 

if ((nargin >= 1) & (isstr(varargin{1})))
    % Callback invocation: 'ButtonDown'
    feval(varargin{:});
    return;
end

if (nargin < 1)
    GETPT_AX = gca;
    GETPT_FIG = get(GETPT_AX, 'Parent');
else
    if (~ishandle(varargin{1}))
        error('First argument is not a valid handle');
    end
    
    switch get(varargin{1}, 'Type')
    case 'figure'
        GETPT_FIG = varargin{1};
        GETPT_AX = get(GETPT_FIG, 'CurrentAxes');
        if (isempty(GETPT_AX))
            GETPT_AX = axes('Parent', GETPT_FIG);
        end

    case 'axes'
        GETPT_AX = varargin{1};
        GETPT_FIG = get(GETPT_AX, 'Parent');

    otherwise
        error('First argument should be a figure or axes handle');

    end
end

% Bring target figure forward
figure(GETPT_FIG);

% Remember initial figure state
buttonDownFcn = get(GETPT_FIG, 'WindowButtonDownFcn');
buttonUpFcn = get(GETPT_FIG, 'WindowButtonUpFcn');
interruptible = get(GETPT_FIG, 'Interruptible');
pointer = get(GETPT_FIG, 'Pointer');

% Remember the ButtonDownFcn's for all axes and images
figaxes = findobj(GETPT_FIG, 'type', 'axes');
figimages = findobj(GETPT_FIG, 'type', 'image');
childButtonDownFcns = get([figaxes;figimages], 'ButtonDownFcn');
set([figaxes;figimages], 'ButtonDownFcn', '');

% Set up initial callbacks for initial stage
set(GETPT_FIG, 'WindowButtonDownFcn', 'getpt(''ButtonDown'');');

% Initialize the lines to be used for the drag
markerSize = 9;
GETPT_H1 = line('XData', [], ...
                  'YData', [], ...
                  'Visible', 'off', ...
                  'Clipping', 'off', ...
                  'Color', 'c', ...
                  'LineStyle', 'none', ...
                  'Marker', '+', ...
                  'MarkerSize', markerSize, ...
                  'EraseMode', 'xor');

% We're ready; wait for the user to do the drag
% Wrap the call to waitfor in eval(try,catch) so we'll
% have a chance to clean up after ourselves.
errCatch = 0;
eval('waitfor(GETPT_H1, ''UserData'', ''Completed'');', 'errCatch=1;');

% After the waitfor, if GETPT_H1 is still valid
% and its UserData is 'Completed', then the user
% completed the drag.  If not, the user interrupted
% the action somehow, perhaps by a Ctrl-C in the
% command window or by closing the figure.

if (errCatch == 1)
    errStatus = 'trap';
    
elseif (~ishandle(GETPT_H1) | ...
            ~strcmp(get(GETPT_H1, 'UserData'), 'Completed'))
    errStatus = 'unknown';
    
else
    errStatus = 'ok';
    x = get(GETPT_H1, 'XData');
    y = get(GETPT_H1, 'YData');
    % If no point was selected, return pt as an empty matrix.
    if (isempty(x) | isempty(y))
       pt = []
    else
       pt = [x,y];
    end
    if (isempty(y))
        y = zeros(0,1);
    end
end

% Delete the animation objects
if (ishandle(GETPT_H1))
    delete(GETPT_H1);
end

% Restore the figure state
if (ishandle(GETPT_FIG))
    set(GETPT_FIG, 'WindowButtonDownFcn', buttonDownFcn, ...
                     'WindowButtonUpFcn', buttonUpFcn, ...
                     'Pointer', pointer, ...
                     'Interruptible', interruptible);
    % Restore the children's button down functions
    set([figaxes;figimages], {'ButtonDownFcn'}, childButtonDownFcns);
end

% Clean up the global workspace
clear global GETPT_FIG GETPT_AX GETPT_H1
clear global GETPT_PT1 

% Depending on the error status, return the answer or generate
% an error message.
switch errStatus
case 'ok'
    % No action needed.
    
case 'trap'
    % An error was trapped during the waitfor
    error('Interruption during mouse point selection.');
    
case 'unknown'
    % User did something to cause the point selection to
    % terminate abnormally.  For example, we would get here
    % if the user closed the figure in the middle of the selection.
    error('Interruption during mouse point selection.');
end


%--------------------------------------------------
% Subfunction ButtonDown
%--------------------------------------------------
function ButtonDown

global GETPT_FIG GETPT_AX GETPT_H1

[x,y] = getcurpt(GETPT_AX);

set(GETPT_H1, ...
        'XData', x, ...
        'YData', y, ...
        'Visible', 'on');
     
set(GETPT_H1, 'UserData', 'Completed');

        
