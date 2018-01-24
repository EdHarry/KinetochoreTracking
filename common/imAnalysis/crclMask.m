function [cMask]=crclMask(imL,imW,centerYX,radius,nRands)
%
% DESCRIPTION: Creates mask of circle with given center and radius
%
% SYNOPSIS: [cMask]=crclMask(im,centerYX,radius)
%
% INPUT: 
%        im         : image matrix
%        centerYX   : (row, col) index for the circle center
%        radius      : circle radius
%
% OUTPUT: 
%        cMask   : im-sized matrix with values from 0-1, depending on
%        whether the pixels fall in the circle (1) or not (0)
%                     
% MATLAB VERSION (originally written on): 7.0.1.24704 (R14) Service Pack 1 Windows_NT
%
% USERNAME: kathomps
% DATE: 2-May-2006
%
% COMMENTS: this function works by sampling each pixel around the circle
% edge with random numbers and calculating how much of its area falls
% within the circle. pixels entirely within the circle are included (but
% not sampled).

cMask=zeros(imL,imW);

% given the center point, get coordinates of all the pixels in the square
% circumscribing the circle
[R C]=meshgrid(centerYX(1)-radius:centerYX(1)+radius, ...
    centerYX(2)-radius:centerYX(2)+radius);

% row and column indices
coordsYX=[R(:) C(:)];

% get array index from xy-coordinates
[indexList]=xy2index(coordsYX(:,2),coordsYX(:,1),imL,imW,1);

if sum(isnan(indexList))~=0
    % if any part of the square falls outside the image boundary,
    % return empty mask and message
    disp('the circle with the following center and radius fell outside the image');
    centerYX
    radius
else
    % if entire circle fits in image, return sum of pixel intensities
    % in sumTotal

    % get yx-coordinates of top left, bottom left, top right, and
    % bottom right corners of each pixel. eg) for pixel 1 with indices
    % (1,1): TL=(0.5,0.5), BL=(1.5,0.5), TR=(0.5,1.5), BR=(1.5,1.5)
    TL=[coordsYX(:,1)-0.5 coordsYX(:,2)-0.5];
    BL=[coordsYX(:,1)+0.5 coordsYX(:,2)-0.5];
    TR=[coordsYX(:,1)-0.5 coordsYX(:,2)+0.5];
    BR=[coordsYX(:,1)+0.5 coordsYX(:,2)+0.5];

    % count the number of corners of each pixel inside the circle
    cornersIn=double(dist2Pt(TL,centerYX(:))<radius)+...
        double(dist2Pt(BL,centerYX(:))<radius)+...
        double(dist2Pt(TR,centerYX(:))<radius)+...
        double(dist2Pt(BR,centerYX(:))<radius);

    % for pixels with 4 corners in, fraction of area in circle is 1.0
    fractionInCrcl=double(cornersIn==4);

    % randPts contains 100 [R C] indices for each pixel
    randPts=rand(sum(cornersIn~=0 & cornersIn~=4),2,nRands)+repmat(coordsYX(cornersIn~=0 & cornersIn~=4,:)-0.5,[1 1 nRands]);

    % get distance of each random point to the center of the circle
    dist=dist2Pt(randPts,centerYX);

    % fraction of random points falling in the circle (corresponds to
    % fraction of the pixel area falling in circle)
    fractionInCrcl(cornersIn~=0 & cornersIn~=4)=sum((dist<=radius),2)/nRands;

    cMask(indexList)=fractionInCrcl;

    % multiply intensities times the fraction of area within circle
    
end

