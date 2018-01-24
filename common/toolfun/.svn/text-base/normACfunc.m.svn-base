function[acnx, acny]=normACfunc(image,verbose)
%
%creates normalized spatial autocorrelation function (of mean zero) of
%image with wraparound, one in x-, one in y-direction. The maximum lag is
%ten times the longer side of the image, or one-half the number of elements
%in the image, whichever is smaller (the acf is symmetric!)
%
% SYNOPSIS: [acnx, acny]=normACfunc(image, verbose)
%
% INPUT: image = 2D-array
%        verbose = optional argument. -1 does not do anything, 0 writes
%                  progress to command line, 1 plots a figure of the
%                  autocorrelation. Default: 1       
%        
% OUTPUT: acnx = normalized autocorrelation function of image projcted in x
%               -direction
%        acny = normalized autocorrelation function of image projcted in y
%               -direction
%last modified by Dinah 07/08/2005

if nargin < 2 || isempty(verbose)
    verbose = 1;
end

[sizex,sizey]=size(image);
dbimage=double(image);
dbimageCon=dbimage';

% number of lags: 10x the longer side of the image, or 0.5x the number of
% pixels, whichever is smaller
nLags = min(10*max(sizex, sizey),ceil(prod([sizex,sizey])/2));

%the mean m ist subtracted from the image further down to ensure that the
%mean of the autocorrelation is zero
m=mean(dbimage(:));

im_xproject=dbimage(:)-m;
im_yproject=dbimageCon(:)-m;

%initialize
shiftedvecx=im_xproject;
shiftedvecy=im_yproject;
[acnx, acny] = deal(zeros(nLags,1));



%variable for progress monitoring
ct=round(nLags/10);

% calculate the normalization parameter already
normx = sum(im_xproject.^2);
normy = sum(im_yproject.^2);

%loop over r = spatial shift in pixels
for r=0:(nLags-1)
    if verbose > -1 && (mod(r,ct)==0) 
        perc = 10*round(r/ct);
        disp(['still computing... progress = ',num2str(perc),'%']);
    end
    %shiftedvec equals original vector shifted by r, with wraparound
    shiftedvecx(1:end-r) = im_xproject((1+r):end);
    shiftedvecx(end-r:end) = im_xproject(1:(1+r));
    
    shiftedvecy(1:end-r) = im_yproject((1+r):end);
    shiftedvecy(end-r:end) = im_yproject(1:(1+r));
    
    %acn =  normalized autocorrelation function int(f(s)*f(s+r))/int(f(s)^2)
    %the total mean of the vector is subtracted to ensure that the
    %autocorrelation function has mean zero (and thus drops to zero for
    %perfectly uncorrelated images)
    acnx(r+1)=sum(im_xproject .* shiftedvecx) / normx;
    acny(r+1)=sum(im_yproject .* shiftedvecy) / normy;

end % of for

if verbose == 1
    account=0:(nLags-1);
    plot(account,acnx,'-r.');
    my=min(min(acnx),min(acny));
    mx=4*max(sizex,sizey);
    axis([-10 mx my 1]);
    hold on
    plot(account,acny,'b-.');
    hold off
end
