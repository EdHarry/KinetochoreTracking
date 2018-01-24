function varargout = histND(varargin)
%HISTND generates an N-dimensional histogram
%
% SYNOPSIS     histND(data)
%              histND(data,nBins)
%              histND(data,'auto',factor)
%        [...] = histND(data, ... , display)
%        [counts] = histND(...)
%        [counts, bins] = histND(...)
%        [counts, bins, spotCoords] = histND(...)
%
% INPUT     data : n-by-d array of n data points in d dimensions
%           nBins (opt): scalar or 1-by-d array of the number of bins
%           'auto' (opt): calculates an optimal bin width that can be
%               modified by multiplication with the optional input argument
%               factor
%           display (opt): if no argument 'display' is given, histND plots
%               into the current axes if no output is requested. With the
%               argument, it will open a figure, or imaris, and plot
%               regardless of the number of output arguments. If the input
%               has too many dimensions for visualization, nothing is
%               displayed.
%               Display options:
%                   'figure' - shows a bar plot for 1D data, or a 2D image
%                               of the counts with imtool (d: 1-2)
%                   'imaris' - shows the counts in Imaris (d: 1-5)
%                   'imarisPlot3' - shows the counts with the values
%                               overlaid in Imaris (d: 1-3)
%
%                   For displaying a 3D bar plot, call 'bar3(counts)'
%
% OUTPUT    counts : number of counts for each bin
%           bins   : position of the center of each bin
%           spotCoords : coordinates of the individual data points relative
%                    to the bins
%               
% c: jonas, 8/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=====================
% TEST INPUT
%=====================

% argument 1
if nargin < 1 || isempty(varargin{1}) || ~isnumeric(varargin{1})
    error('histND requires a non-empty numeric data array!')
else
    data = varargin{1};
    if ndims(data) > 2
        error('only 2D-arrays can be processed by histND!')
    end

    [nPoints, nDims] = size(data);
    minData = min(data,[],1);
    [maxData] = max(data,[],1);
    % loop through dimensions to check if there have been several maxIdx.
    % Check individually for each col of data, because there are several
    % different maxima
    maxIdx = [];
    for d=1:nDims
        maxIdx = [maxIdx; find(data(:,d) == maxData(d)) + nPoints * (d-1)];
    end
end

% argument 2
if nargin < 2 || isempty(varargin{2})
    % default: 10 bins along each direction
    nBins = 10 * ones(1,nDims);
else
    % decide whether nBins or 'auto'
    if isnumeric(varargin{2})
        nBins = varargin{2};
        nBins = nBins(:)';

        % check dimension
        switch length(nBins)
            case 1
                nBins = nBins * ones(1,nDims);
            case nDims
                % all is well
            otherwise
                error(['nBins has to be either a scalar or a 1-by-d'...
                    'array, where d is the number of cols of data!'])
        end

    elseif ischar(varargin{2}) & strcmpi(varargin{2},'auto')
        % prepare nBins for automatic calculation
        nBins = -1;

        % check for varargin{3} - is it factor
        if nargin > 2 && ~isempty(varargin{3}) && ~ischar(varargin{3})
            factor = varargin{3};
            factor = factor(:)';
            % remove factor from input argument list
            varargin(3) = [];

            % check dimension
            switch length(factor)
                case 1
                    factor = factor * ones(1,nDims);
                case nDims
                    % all is well
                otherwise
                    error(['factor has to be either a scalar or a 1-by-d'...
                        'array, where d is the number of cols of data!'])
            end
        else
            factor = ones(1,nDims);
        end

    else % bizarre input
        error('histND requires an array of bins or ''auto'' as second input argument!')
    end
end

% set display option - we can't use nargin anymore, because we might
% have removed factor from input argument list
% 1-d display bar plot (1)
% 2-d display with imtool(M7)/imshow(M6)/imaris (2/2/3)
% 3-d, 4-d display with imaris with array and points (4)
% 3-d, 4-d, 5-d display with imaris as array only (3)
% 6-d cannot display

% default: no display
display = 0;

if length(varargin) < 3 || isempty(varargin{3})
    % if no display requested, set only if no output arguments
    if nargout == 0 && nDims < 3
        % we can display. Set option later
    elseif nargout > 0
        % we don't care
    else
        % exit here, because there is no display and no output
        warning('too many dimensions. Cannot display output')
        return
    end
elseif ischar(varargin{3})
    switch varargin{3}
        case 'figure'
            if nDims < 3
                display = nDims;
            else
                % warn and return if no output
                warning('too many dimensions. Cannot display output')
                if nargout == 0
                    return
                end
            end
        case 'imaris'
            if nDims < 6
                % if nDims = 1, imaris might not be the best choice, but it
                % is one the user made
                display = 3;
            else
                % warn and return if no output
                warning('too many dimensions. Cannot display output')
                if nargout == 0
                    return
                end
            end
        case 'imarisPlot3'
            if nDims < 4
                % if nDims = 1, imaris might not be the best choice, but it
                % is one the user made
                % We can't have more than 3 dims, because time cannot be a
                % float!!
                display = 4;
            else
                % warn and return if no output
                warning('too many dimensions. Cannot display output')
                if nargout == 0
                    return
                end
            end
        otherwise
            error('bad option for display mode')
    end
else
    error('bad option for display mode')
end

% check for too many output arguments
if nargout > 3
    error('too many output arguments!')
end

%==========================================================================


%========================
% HISTOGRAM PARAMETERS
%========================

% find nBins if we haven't already
if nBins == -1
    if nPoints < 20
        warning('Less than 20 data points!')
        nBins = ceil(nPoints/4) * ones(1,nDims);
    else


        % create bins with the optimal bin width
        % W = 2*(IQD)*N^(-1/3)
        interQuartileDist = diff(prctile(data,[25,75]));
        binLength = 2*interQuartileDist*length(data).^(-1/3).*factor;

        % number of bins: divide data range by binLength. Err towards
        % larger bins, but make sure that they contain all the data!
        nBins = floor((maxData-minData)./binLength);
        binLength = (maxData - minData)./nBins;
        
        % nBins can be NaN for discrete distributions. Take care of that
        % here
        for d = 1:nDims
            if ~isfinite(nBins(d))
                nBins(d) = length(unique(data(:,d)));
            end
        end

    end
else
    % calculate binLength
    binLength = (maxData - minData)./nBins;
end

%========================



%========================
% HISTOGRAM
%========================

% calculate which bins the data points fall into. SubPixPos will be the
% position of the data point in the bin-coordinates (for plotting spots in
% imaris)
w = warning;
warning('off','RETURNRIGHTVECTOR:SQUAREINPUT');
[subPixPos, whichBin] = wraparound(data,[minData;minData+binLength]);
warning(w);
% Normally, if a value would fall right between two bins, it will be
% counted for the bin on the right. However, this should not be the case
% for the maximum that should belong to the last bin
whichBin(maxIdx) = whichBin(maxIdx) - 1;
% stretch subPixPos so that it lies between 0 and 1
subPixPos = (subPixPos-repmat(minData,nPoints,1))./repmat(binLength,nPoints,1);
% transform to image zero
subPixPos = subPixPos + whichBin + 0.5;
% first bin is 1,1,1
whichBin = whichBin + 1;


% for counts: countEntries of whichBin - which is nothing else than an
% array of pixel coordinates
c = mat2cell(whichBin,nPoints,ones(1,nDims));
countIdx = sub2ind(nBins,c{:});
countIdx(isnan(countIdx)) = [];
[countIdx,n] = countEntries(countIdx);
if nDims == 1
    counts = zeros(nBins,1);
else
    counts = zeros(nBins);
end
counts(countIdx) = n;

%======================


%======================
% DISPLAY
%======================

switch display
    case 0
        if nargout == 0
            % display. We know from above that nDims < 3
            switch nDims
                case 1
                    bar(gca,bins,counts,'hist');
                case 2
                    imshow(counts);
            end
        end
    case 1
        % bar plot
        figure('Name','Histogram')
        bar(gca,bins,counts,'hist')

    case 2
        % imshow or imtool
        if strmatch('6',version)
            uiViewPanel;
            imshow(counts,[]);
        else
            try
                imtool(counts,[]);
            catch
                figure,
                imshow(counts,[]);
            end
        end

    case 3
        % imarisShowArray
        imarisShowArray(counts)

    case 4
        % imarisPlot3
        plotData.XYZ = subPixPos;
        imarisPlot3(plotData,[],counts)
end


%======================
% ASSIGN OUTPUT
%======================


switch nargout
    case 0
        varargout = [];
    case 1
        varargout{1} = counts;
    case 2
        varargout{1} = counts;
        for d = nDims:-1:1
            bins{d} = ...
                [minData(d)+0.5*binLength(d):binLength(d):maxData(d)]';
        end
        varargout{2} = bins;
    case 3
        varargout{1} = counts;
        for d = nDims:-1:1
            bins{d} = ...
                [minData(d)+0.5*binLength(d):binLength(d):maxData(d)]';
        end
        varargout{2} = bins;
        varargout{3} = subPixPos;
end