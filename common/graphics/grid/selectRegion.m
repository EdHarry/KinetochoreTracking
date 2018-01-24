function [Gm,Wm,G,W,bw,bwMask,mnWindows,nWindows,iBwMask]=selectRegion(img,GridSize,WindowSize,iBwMask)
% SELECTREGION  allows the user to define a grid of subwindows and to select some of them by drawing a region of interest onto the image
%
% SYNOPSIS      [Gm,Wm,G,W,bw,bwMask,mnWindows,nWindows]=selectRegion(img,GridSize,WindowSize)
%
% INPUT         img        : image from which a region of interest has to be chosen
%               GridSize   : size of the grid of sub-window centers
%               WindowSize : size of the subwindows centered on the grid defined by GridSize
%               iBwMask    : if a B&W mask is passed, the user does not need to draw a region of interest to
%                            select a subset of the grid; instead, the mask is used.
%
%               Note: selectRegion uses the function FRAMEWORK to generate the sub-window grid. 
%
% OUTPUT        Gm         : centers of all sub-windows fully contained in the region of interest
%               Wm         : vertex coordinates of all subwindows fully contained in the region of interest
%               G          : centers of all sub-windows in the image independent of selection
%               W          : vertex coordinates of all subwindows independent of selection
%               bwMask     : boolean matrix defining which subwindows are in the ROI and which ones are outside
%               mnWindows  : number of subwindows in y- and x-direction independent of selection
%               nWindows   : number of windows contained in the region of interest
%               iBwMask    : draen (or input) B&W mask
%
% DEPENDENCES   selectRegion uses { framework }
%               selectRegion is used by { }
%
% Aaron Ponti, 2001

if nargin==3
    iBwMask=[];
end

% Check iBwMask
if ~isempty(iBwMask)
    if size(iBwMask)~=size(img) | size(iBwMask,3)~=1 | (length(find(iBwMask~=0))+length(find(iBwMask~=1))~=prod(size(img)))
        error('iBwMask is not valid');                
    end
end

% Initialize grids
[G,W,S,mnWindows]=framework(size(img),GridSize,WindowSize);
clear S;

% Display current image
figure;
imshow(img,[]);
set(gcf,'NumberTitle','off');
hold on;
    
% Display grid centers
plot(G(:,2),G(:,1),'bo');

if isempty(iBwMask)
    
    % Title
    set(gcf,'Name','Please draw region of interest');

    % Select region of interest
    try
        iBwMask=roipoly; 
        iBwMask=double(iBwMask);
    catch
        Gm=[];Wm=[];G=[];W=[];bw=[];
        bwMask=[];mnWindows=[];nWindows=[];iBwMask=[];
        close(gcf);
        disp('Aborted.');
        return
    end

else
    
    % Title
    set(gcf,'Name','Grid');
    
end

% Erase windows not contained in the polygon
result=(W(1,3)-W(1,1)+1)*(W(1,4)-W(1,2)+1);
counter=0;

toBeDel=[];
for j=1:size(W,1)
   testArea=iBwMask(W(j,1):W(j,3),W(j,2):W(j,4));
   if sum(testArea(:))~=result
      % The window is NOT completely inside the selected area
      counter=counter+1;
      toBeDel(counter)=j;
   end
end

% Delete all sub-windows outside the polygon
Gm=G;
Wm=W;
if ~isempty(toBeDel)
   Gm(toBeDel,:)=[];
   Wm(toBeDel,:)=[];
end

% Number of actual windows
nWindows=size(Wm,1);

% Create black and white image
bw=zeros(size(img));
for i=1:size(Wm,1)
   bw(Wm(i,1):Wm(i,3),Wm(i,2):Wm(i,4))=1;
end

% Create mask for selection of the subwindows
bwMask=ones(mnWindows(2),mnWindows(1));
bwMask(toBeDel)=0;
bwMask=bwMask';

% Show result
imshow(bw,[]);
plot(Gm(:,2),Gm(:,1),'bo');
pause(1)
close;
