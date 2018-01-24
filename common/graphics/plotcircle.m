function lineHandle = plotcircle(varargin)
%PLOTCIRCLE draws a circle
%
% SYNOPSIS: lineHandle = plotcircle(center,radius,lineSpec)
%           lineHandle = plotcircle(center,radius,pn,pv)
%           lineHandle = plotcircle(ah,...)
%
% INPUT center (optional): n-by-2 list of centers. Default: [0,0]
%       radius (optional): n-by-1 list of radii. Default: [1]
%       lineSpec (optional): see 'doc linespec'
%       pn,pv (optional): line property name/property value pair(s). See 
%               'doc line_props', or 'doc set' for details
%
% OUTPUT lineHandle: list of handles
%
% REMARKS
%
% created with MATLAB ver.: 7.5.0.342 (R2007b) on Windows_NT
%
% created by: Thai-Hang Nguyen
% DATE: 28-Mar-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=================
%% TEST INPUT
%=================

%assign defaults
def_center = [0,0];
def_radius = [1];
nPoints = 500; % number of points on a circle


% check for axes handle
if nargin > 0 && ~isempty(varargin{1}) && length(varargin{1})==1 && ishandle(varargin{1}) && strmatch(get(varargin{1},'Type'),'axes')
    axesHandle = varargin{1};
    varargin(1) = [];
else
    % default: take current axes
    axesHandle = gca;
end

numArgIn = length(varargin); % assign the number of input elements to numArgIn

%check for center
if numArgIn < 1 || isempty(varargin{1})
    % assign default (expand to match number of radii)
    center = def_center;
else
    % assign center
    center = varargin{1};
end % end of if statement:check center

% check for radius
if numArgIn < 2 || isempty(varargin{2})
    %assign default
    radius = def_radius;
else
    if any(varargin{2} < 0)
        error ('drawCircle needs a positive radius');
    else
        %assign radius
        radius = varargin{2};
    end

end % end of if statement:check radius

%number of radii
radius = radius (:);
nRadii = length (radius);
%number of centers
[nCenters,nDims] = size(center);

% expand centers and radii if necessary
if nCenters == 1 && nRadii > 1
    center = repmat(center,nRadii,1);
    nCenters = nRadii;
end
if nCenters > 1 && nRadii == 1
    radius = repmat(radius,nCenters,1);
    nRadii = nCenters;
end
if nCenters ~= nRadii || nDims ~= 2
    error('drawCircle requires center to have the same number of rows as radius and two columns')
end




% define default options (none)
plotOpt = '';
setCell = '';
switch max(numArgIn -2,0)
    case 0
        %no line options entered, keep default
    case 1
        %check if a numerical input and throw an error
        if isnumeric(varargin{3})
            error('String or cellarray types are accepted.')
        end
        % check if a string input
        if ischar(varargin{3})
            plotOpt = varargin{3};
        end


    otherwise
        if isEven(numArgIn-2)
            %  if 2 cellarrays, do additional checking
            if iscell(varargin{3}) && iscell(varargin{4})
                %check for consistency of the properties inputs when 2 pn
                %should also have 2 corresponding pv
                [pnrow,pncol] = size(varargin{3});
                [pvrow,pvcol] = size(varargin{4});
                if pncol ~= pvcol
                    error('missing a pn or a pv')
                end

                if pvrow ~= nRadii
                    error('The number of line options has to be equal to the number of radii (number of circles to be drawn)')
                end
            end % if there are two cellArrays

            setCell = varargin(3:end);
            %disp(setcell);

        else
            error('a pn or a pv input is missing')
        end

end


%==================
%% CALCULATE CIRCLE
%==================
% Take different points following the equations and connect them to create
% a circle
% x = Rcos(theta) + x(centre)
% y = Rsin(theta) + y(centre)
% theta is from 0 to 2pi

%theta = 0:(2*pi)/(nPoints-1):2*pi;
theta = linspace(0,2*pi,nPoints); % takes the number of points around the circle

% * multiply arrays //  .* is multiply elements one by one
% Here multiply 2 1D arrays to get a 2D array where one row is the list of
% x or y coordinates of one circle

x = radius*cos(theta) + repmat(center (:,1), 1, nPoints);  % takes the first column and expand nPoints times
y = radius*sin(theta) + repmat(center (:,2), 1, nPoints);




%==================
%% DRAW CIRCLE
%==================

% plot all at once
lh = plot (axesHandle, x',y',plotOpt);

% set additional properties only if there are any
if ~isempty(setCell)
    set(lh, setCell{:});
end
axis(axesHandle,'equal')

% only return line handle if requested
if nargout > 0
    lineHandle = lh;
end



