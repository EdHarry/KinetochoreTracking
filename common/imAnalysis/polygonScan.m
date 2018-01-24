function [polyScanOut,polygonIn] = polygonScan(imageIn,polygonIn,showPlots)

%[polyScanOut,polygonIn] = polygonScan(imageIn,polygonIn,showPlots)
%
%Performs a "polygon/line scan" of the input image. This returns the image
%values at pixels which are traversed by the input polygon.
%If no polygon is input, the user is asked to click on the image to create
%one.
%Warning: If the polygon intersects itself, weird shit is going to happen.
% 
% Input:
% 
%   imageIn - the image you want to perform a linescan on.
% 
%   polygonIn   -   Optional. The polygon / line describing the scan to
%                          perform. If left empty, the user enters one
%                          interactively.
% 
%   showPlots   -   If true, a plot of the resulting scan is displayed.
%
%
%
%   Output:
% 
%   polyScanOut     -   This is a Mx2 array containing the values of each
%                                 pixel along the line (2nd column), and
%                                 its relative position along the length of
%                                 each side of the polygon (1st column).
%                                 Points on side 1 will go from 1.000 to
%                                 1.999, points on side to will go from
%                                 2.000 to 2.999 and so on. This means if
%                                 the different sides are different lengths
%                                 you need to be carefu comparing theml!
% 
% 
% 
%  polygonIn        -    This is the 2xM vector of the vertex positions of
%                               the polygon used to scan the image.
% 
% 
% 
% 
% 
% 
 %Hunter Elliott, 3/2009


if nargin < 1 || isempty(imageIn)
    error('Must input an image!')
end



%If no polygon input, allow the user to create one interactively.
if nargin < 2 || isempty(polygonIn)
    
    figHan = figure;        
    imagesc(imageIn);
    axis image,colormap gray
    polyHan = impoly(get(figHan,'CurrentAxes'),'Closed',0);    
    polygonIn = getPosition(polyHan)';
    
end

if nargin < 3 || isempty(showPlots)
    showPlots = false;
end



%Get the pixels this polygon traverses.
polyMat = drawpoly(polygonIn,size(imageIn));

%Get the position along the line
linePos = polyMat(polyMat>0);
[linePos,sortInd] = sort(linePos);
%Get the image values at these points
imVals = cast(imageIn(polyMat>0),'double');
imVals = imVals(sortInd);
polyScanOut = [linePos,imVals];

if showPlots
    figure
    plot(polyScanOut(:,1),polyScanOut(:,2),'.-')
    xlabel('Position along polygon side(s)');
    ylabel('Intensity')
end