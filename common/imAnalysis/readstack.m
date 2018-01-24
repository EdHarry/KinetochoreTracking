function readstack(img_first,img_last)

% readstack loads a metamorph stk into a structure in the workspace using metaTiffRead
%
% SYNOPSIS   readstack
%
% INPUT      metamorph stack
%
% OUTPUT     none     :    
% 
% REMARKS    
%       
%
%
% DEPENDENCES   readstack uses {metaTiffRead, readstackSlider}
%               readstack is used by {nothing}
%
% Johan de Rooij, June 10th, 2005

% get the desired stack read in matlab as a structure, determine length at
% the same time.
if (nargin == 0)
    img_first = 1;
    img_last = 10000;
end

if (nargin==1)  img_last = img_first;  end


[Stack]=metaTiffReadJR(img_first,img_last);
maxnr = length(Stack);

% Draw a new figure on the screen, create a handle to it.
hFigure = figure;

% Draw the frame counter in the figure; it is identified by the tag picturecount
FrameCounterHandle = uicontrol ('Style', 'text',...
                                'Units', 'normalized',...
                                'Tag', 'framecount',...
                                'Position', [0.5,0.93,0.05,0.06]);

% Set the frame counter to the first image number
set (FrameCounterHandle, 'String',img_first);

% Draw the slider in the figure; it is identified by the tag pictureslide and calls
% the function readstackSlider when moved
sliderHandle = uicontrol ('Style', 'slider', ...
                          'Units', 'normalized', ... 
                          'Value', 1/(maxnr), ...
                          'Min', 1/(maxnr), ...
                          'Max', 1, ...
                          'SliderStep', [1/maxnr 5/maxnr], ...
                          'Callback', 'readstackSlider', ...
                          'Tag', 'pictureslide', ...
                          'Position', [0.02,0.02,0.05,0.9]);



% Show the first frame on the screen in the current figure
frame1 = double(Stack(1).data);
hold on;
%imshow (image, []), title (num2str (firstImage));
imshow (frame1, [min(frame1(:)) max(median(frame1(:)))*2]);
hold off;

%for use in callback..
handles.maxnr = maxnr;
handles.Stack = Stack;
handles.img_first = img_first;
guidata(hFigure,handles);

% et voila! (with thanks to Andre Kerstens and Colin Glass)
