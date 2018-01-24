function varargout = arbitraryGrid(e_1, e_2, e_3, zero, v1, v2, v3, imSize)
% ARBITRARYGRID calculates an arbitrarily oriented, regularly spaced grid
% in 3D (for 2D: input [0,0,1] for e_3) for reading values from an image with
% an interpolator
%
% SYNOPSIS varargout = arbitraryGrid(e_1, e_2, e_3, zero, v1, v2, v3, imSize)
%         
% INTPUT   e_1,e_2,e_3 : unit vectors of the grid. Spacing of the grid
%                        depends on the length of the vectors
%          zero        : origin of the grid 
%          v1,v2,v3    : for each argument: 
%                        - specify numeric array [min max]
%                          to indicate the extent of the grid. [-2 3] will
%                          generate 6 grid points in the given direction,
%                          with the zero-point as the third element
%                        - specify inf for either min or max or both if
%                          the extent should be chosen to 
%                          include all of the original data. This option
%                          requires the input of the size of the image,
%                          imSize
%                        **** later: option 'tight' , which needs imSize
%
% OUTPUT   listOfCoords OR [xGrid,yGrid,zGrid]. If only one output argument
%                       is asked for, the grid will be given as a list of
%                       coordinates. If several output arguments are asked
%                       for, the output is given in ndgrid-style
%                       (required for interpolation)
%
%          
% c: jonas 11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================
% PARSE INPUT
%================

% number of input arguments
error(nargchk(7, 8, nargin, 'struct'));

% make sure unit vectors are row vectors
e_1 = returnRightVector(e_1, 1, 'r');
e_2 = returnRightVector(e_2, 1, 'r');
e_3 = returnRightVector(e_3, 1, 'r');
zero = returnRightVector(zero, 1, 'r');
v1 = returnRightVector(v1, 1, 'r');
v2 = returnRightVector(v2, 1, 'r');
v3 = returnRightVector(v3, 1, 'r');

%================


%=========================
% PREPARE GRID-GENERATION
%=========================

%calculate grid region here
allExtent = [v1;v2;v3];
infExtent = abs(allExtent) == inf;

if any(infExtent(:))
    % calculate where the corners of the original image would be within the
    % new grid
    if nargin < 8 || isempty(imSize) 
        error('If you want ARBITRARYGRID to calculate the extent of the grid, you have to specify the image size ''imSize''');
    end
    e_matrix = [e_1',e_2',e_3'];
    [u,v,w] = ndgrid([1,imSize(1)],[1,imSize(2)],[1,imSize(3)]);
    imageCorners = [u(:)';v(:)';w(:)']-repmat(zero',1,8);
    requiredSize = ceil(e_matrix\imageCorners);
    
    % find minimum and maximum extent for the three new dimensions
    minSize = min(requiredSize,[],2);
    maxSize = max(requiredSize,[],2);
    
    % fill in allExtent, then separate into v1-3 again
    reqSize = [minSize,maxSize];
    allExtent(infExtent) = reqSize(infExtent);
    v1 = allExtent(1,:);
    v2 = allExtent(2,:);
    v3 = allExtent(3,:);
end

% count number of elements along each dimension. Don't forget to add 1 to
% include the first element
num1 = floor(diff(v1))+1;
num2 = floor(diff(v2))+1;
num3 = floor(diff(v3))+1;

% preassign output array
listOfCoords = zeros(num1*num2*num3,3);

%========================


%========================
% GENERATE GRID
%========================

% loop through dimensions so that list of coordinates is ordered according
% to their indices
k = 1;
for z = v3(1):v3(2)
    for y = v2(1):v2(2)
        for x = v1(1):v1(2)
            
            % linear combination of unit vectors
            listOfCoords(k,:) = x*e_1 + y*e_2 + z*e_3 + zero;
            
            % up index by one
            k = k+1;
            
        end
    end
end

%=======================


%=======================
% ASSIGN OUTPUT
%=======================

% return 
if nargout == 1
   varargout{1} = listOfCoords;
end
if nargout > 1
    varargout{1} = reshape(listOfCoords(:,1),num1,num2,num3);
    varargout{2} = reshape(listOfCoords(:,2),num1,num2,num3);
end
if nargout > 2
    varargout{3} = reshape(listOfCoords(:,3),num1,num2,num3);
end