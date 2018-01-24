% MEXFGINTERFACE is the framegrabber to matlab interface
%
% SYNOPSIS [s,pics] = mexFgInterface(command,p1,p2,p3,..)
%
% INPUT: command string , parameters
%        'live'         , p1 = filename string (optional)
%        'open'         , p1 = framegrabber value (0 = PicPortMono)
%        'grab'         , p1 = roi=[top,left,width,height],
%                         p2 = nr of images per camera,
%                         p3 = shutter time in micro seconds
%        'close'
%
% OUTPUT: s    : status value (OK=1, ERROR != 1)
%         pics : matrix [m,n,a];
%                a = picture index
%                matrix [m,n] = (256 colors/grey values) picture
% Note:
% 1. All grab params are optional, if no params are given the last
%    configuration remains active.
% 2. The default values are p1 = full image, p2 = one image per camera
%    p3 = maximal shutter time
% 3. Without params the grab command is faster.
%
% 4. 'live' start the live image application. Use the following commands when started:
%    left mousebutton = start live
%    right mousebutton = stop live
%    + keypad = increase shuttertime
%    - keypad = decrease shuttertime
%    space = switch cameras
%    s = save current image to disk