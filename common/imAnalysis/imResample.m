function  dataOut = imResample(dataIn, stepIn, stepOut)
%IMRESAMPLE resamples a n-D array
%  The input data is assumed to be uniformely spaced in any of the
%  dimensions, but the size can vary between dimensions. If the stepsizes
%  do not allow for an array of the same 'size' (e.g. in um for an image),
%  the output image will be smaller, and aligned at [0,0] in matlab
%  coordinates. Transpose the image first, if you measured in image coords.
%  Be careful: Resampling changes the zero position of your image (it lies
%  at -0.5 pixels, and this changes with resampling)
%  If the ratio between the pixelsizes is within 1e-14 of an integer, the
%  ratio will be treated as integer.
%
%  SYNOPSIS dataOut = imResample(dataIn, stepIn, stepOut)
%
%  INPUT    dataIn : n-dim array to be resampled
%           stepIn : pixelsize of the input data (length: n)
%           stepOut: pixelsize of the output data (length: n)
%               if stepOut is larger than stepIn, there will be
%               downsampling.
%
%  OUTPUT   dataOut: n-dim resampled array
%
%  REMARKS  currently, the code will crash with integer ratios and
%           non-commensurate sizes. Solution: append the input with NaNs,
%           and use nansum. Then, multiply the border pixel line so that
%           there is no obvious resampling artifact
%
% c: 10/04 jonas
%    02/06 kathryn - changed mean to sum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===================
% TEST INPUT
%===================

% check for correct number of input arguments
msg = nargchk(3, 3, nargin, 'struct');
if ~isempty(msg)
    disp('wrong number of input arguments')
    help imResample
    return
end

% check length of stepIn, stepOut
oldSize = size(dataIn);
nDims = length(oldSize);

if length(stepIn) ~= nDims
    error('length(stepIn) has to equal length(size(dataIn))');
end

if length(stepOut) ~= nDims
    error('length(stepOut) has to equal length(size(dataIn))');
end

%==================



%=========================================
% LOOP THROUGH DIMENSIONS TO DOWNSAMPLE
%=========================================

% There are two ways for downsampling: In case the ratio of stepsizes is an
% integer, the matrix is reshaped and the sum is taken. In case the ratio
% is a rational number, the cumulative sum is calculated along the
% dimension of interest, then the intensities are read via linear
% interpolation, and finally, the intensity of the new pixel is found by
% taking the difference of the resampled cumulative sum (idea GD).

% calculate pixelRatio. Take care of rounding errors
pixelRatio = stepOut ./ stepIn;
integerRatio = isApproxEqual(round(pixelRatio),pixelRatio,1e-14,'absolute');
pixelRatio(integerRatio) = round(pixelRatio(integerRatio));

% calculate size of the new image. Round down if necessary
newSize = floor(oldSize ./ pixelRatio);

% find diference between old and new size
deltaSize = oldSize - newSize .* pixelRatio;

% init loop
currSize = oldSize;

% now loop through dimensions. Use cumsum approach only when needed (slow
% interp). Shift dimensions because interp1 only works along the first
% dimension.

% delta is the amount of offset from [0,0,0]. If deltaSize/2, the new image
% is centered on the old image. If 0, old and new images are aligned at 0.
delta = 0; % deltaSize/2;

% if this line poses problems, we have to rename the input variable dataOut
% (or just image)
dataOut = dataIn;
clear dataIn

for dim = 1:nDims
    
    % set new size
    currSize(1) = newSize(dim);
    
    % switch between cumsum and reshape
    if integerRatio(dim)
        % integer ratios allow reshape-approach
       
        if pixelRatio(dim) == 1
            % do nothing
        else
        % reshape, so that sum can be taken along dimension 1
        dataOut = reshape(dataOut,pixelRatio(dim),[]);
        
        dataOut = sum(dataOut);
        
        % shape back (size has been updated already)
        dataOut = reshape(dataOut, currSize);
        end
        
    else
        % non-integer ratios require interpolation
        % (this includes up-sampling. In principle, it is possible to
        % upsample in a simple way for 1/integer ratios, but I do not
        % assume that upsampling will be used a lot)
        
        % cumsum
        dataOut = [zeros([1, currSize(2:end)]); cumsum(dataOut,1)];
        
        % interpolate
        dataOut = interp1([0:oldSize(dim)], dataOut, ...
            [0:pixelRatio(dim):oldSize(dim)], 'linear');
        
        % take the difference, adjust intensity for changed pixelsize
        dataOut = diff(dataOut);
        
    end
    
    % shift image by one dimension
    dataOut = shiftdim(dataOut, 1);
    currSize = [currSize(2:end), currSize(1)];
    
end % for dim=1:nDims
