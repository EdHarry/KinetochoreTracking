function readstackSlider
% readstackSlider finds the right image from the current open metamorph
% stack, displays it and puts the frame number above the image. it allows
% the slider to work.
%
% SYNOPSIS       readstackSlider
%
% INPUT          none (it gets values from the slider created in readstack)
%
% OUTPUT         none 
%
% DEPENDENCIES   readstackSlider uses {nothing}
%                                  
%                readstackSlider is used by { readstack }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% johan de rooij        061005          my first, allmost indepently
%                                       created code....


% this is the callback of the slider, created in readstack
% What we do here is:
% - find out which frame the user currently wants to look at
% - show this frame in the figure (created in readstack)

% Look for objects needed for information
% get figure handle
hFigure = gcbo;
handles = guidata(hFigure);

% Get values from the handles
maxnr = handles.maxnr;
Stack = handles.Stack;
img_first = handles.img_first;

sliderHandle = findall (0, 'Tag', 'pictureslide');
FrameCounterHandle = findall (0, 'Style', 'text', 'Tag', 'framecount');

% Get the current value of the slider, so that we know which frame the user wants to process
sliderValue = get (sliderHandle, 'Value');

% calculate framenr
framenr = round (sliderValue * maxnr);

% Write the current frame number in the little window above the slider
set (FrameCounterHandle, 'String', img_first-1+framenr);

% Read the image frame from the Stack structure
CurrentImage = double(Stack(framenr).data);

% Show the frame on the screen in the current figure
hold on;
%imshow (image, []), title (num2str (imageNumber));
imshow (CurrentImage, [min(CurrentImage(:)) max(median(CurrentImage(:)))*2]);
hold off;

% That's it: wait for the next user action
