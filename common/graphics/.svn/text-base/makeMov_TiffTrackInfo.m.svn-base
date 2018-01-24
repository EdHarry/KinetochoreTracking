function []=makeMov_TiffTrackInfo(movie,trackinfo);
% makeMov_TiffTrackInfo makes a QT movie from a matrix of images
% and the trackinfo matrix from Khuloud's tracker
% INPUT:    movie = matrix of images, see below
%           trackinfo = trackInfo matrix from Khuloud's tracker
%
% before running this function, load the TIFFs into a matrix, e.g. by
% running something along these lines:
%
% MovieName = zeros(200,200,300);
% for i=1:301
%   if i<10
%       filename = ['filename00',num2str(i),'.tif'];
%   elseif i<100
%       filename = ['filename0',num2str(i),'.tif'];
%   else
%       filename = ['filename',num2str(i),'.tif'];
%   end
%   imfile = imread(filename,'tif');
%   MovieName(:,:,i) = imfile(1:200,1:200);
% end
%
% The resulting image matrix MovieName has the size [ix,iy,numf], where
% ix = imagesize in x, iy = imagesize in y, numf = number of images
%
% NOTE: if you want to get a movie, then run the code
% MakeQTMovie start 'nameOfMovie'
% before running this function, and afterwards run
% MakeQTMovie finish
% AND uncomment the line 55 below if necessary

[stx,sty] = size(trackinfo);
[ix,iy,numf] = size(movie);

% extract x-positions of tracks
trackinfoX = trackinfo(:,1:8:sty);
% extract y-positions of tracks
trackinfoY = trackinfo(:,2:8:sty);

% display detected and tracked objects with image
for t=1:numf-1
    imshow(movie(:,:,t),[]); hold on;
    % x and y positions
    xpos = nonzeros(trackinfoX(:,t));
    ypos = nonzeros(trackinfoY(:,t));
    if length(xpos)~=length(ypos)
        error('number of x and y positions doesn''t match');
    end
    plot(xpos,ypos,'mo','LineWidth',2);
    pause(0.1);
    
    % if the series should be written into a movie, run the code
    % MakeQTMovie start filename
    % before this function, and uncomment the following line
     MakeQTMovie addaxes
    % then run after this function:
    % MakeQTMovie finish
    
end % of for

end % of function