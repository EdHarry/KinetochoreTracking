function h=vectorFieldPlot(F,handle,polygon,scale)
% vectorFieldPlot scales and displays a vector field or its divergence
%
% SYNOPSIS   h=vectorFieldPlot(F,handle,polygon,scale)
%
% INPUT      F      : either a...
%                      VECTOR FIELD M stored in a (nx4)-matrix of the 
%                         form [y0 x0 y x]n, (where (y0,x0) is the base and (y,x) 
%                         is the tip of the vector), or a...
%                      DIVERGENCE divM stored in a (nx3)-matrix of the form
%                         [y0 x0 div]n, (where (y0,x0) is the base of the vector
%                         and div its divergence).
%            handle : handle of a previous vector field plot for overlapping
%                     (set handle to 0 is you want to draw on a new figure).
%                     This parameter is ignored for divergence field plots.
%            polygon: (optional - pass polygon=[] to disable). The vector
%                     field to be plotted can be cropped to remove vectors 
%                     outside a given region of interest.
%                     To create the polygon use the functions ROIPOLY or
%                     GETLINE. These functions return the polygon vertices
%                     stored in two vectors y and x. 
%                     Set polygon=[x y] to use with vectorFieldPlot.
%            scale  : (optional) defines the scaling factor for F
%
% OUTPUT     h     : handle of the figure
%
% REMARK     Up to 5 vector fields can bo plotted on top of each other (by passing the handle 'h' returned by
%            the previous call to h=vectorFieldPlot(F,h,...) as input) and will be colored in the following order:
%            black,red,blue,magenta,cyan.
%            If they are overlaid onto a figure (e.g. FSM image), the order will be the following:
%            yellow,red,blue,magenta,cyan.
%            Additional vector fields will be blue.
%
% Aaron Ponti 11/28/2002

% If no scale factor is specified, set it to 1
if nargin==2 
    scale=1;
end

% Set all vectors outside the passed polygon to 0
if ~isempty(polygon)
    index=inpolygon(F(:,1),F(:,2),polygon(:,2),polygon(:,1));
    F(find(~index),3)=F(find(~index),1);
    F(find(~index),4)=F(find(~index),2);
end

if size(F,2)==4
    
    % Remove all incomplete lines
    F=F(find(F(:,1)~=0 & F(:,3)~=0),:);
    
    % Use manual scaling instead of quiver scaling which is field-dependent
    if scale~=1
        F(:,3:4)=[F(:,1:2)+scale*(F(:,3:4)-F(:,1:2))];
    end
    
    if handle~=0
        h=figure(handle);
        hold on;
    else
        h=figure;
    end
    
    % Check for MATLAB version - quiver 7.0 is no longer compatible
    v=ver('MATLAB');
    pointPos=findstr(v.Version,'.');
    if ~isempty(pointPos)
        v.Version=v.Version(1:pointPos(1)+1);
    end

    % Plot scaled vector field and change the axis orientation to ij (see help on axis)
    if str2num(v.Version)<7
        quiver(F(:,2),F(:,1),F(:,4)-F(:,2),F(:,3)-F(:,1),0);
    else
        quiver('v6',F(:,2),F(:,1),F(:,4)-F(:,2),F(:,3)-F(:,1),0);
    end
    axis ij
    xlabel('x');
    ylabel('y');
    
    % Get handles of all plots present in the figure
    plotHandles=findall(h,'Type','Line');
    nPlots=length(plotHandles);
    
    % Check is an image exists in the current axis
    imgHandle=findall(h,'Type','Image');
    
    if isempty(imgHandle)
        % Prepare up to five different colors for multiple overlapping
        % plots (this color table is best suited for simple vector field
        % plots; plotted on white background)
        colorTable={'black','red','blue','magenta','cyan'};
    else
        % Prepare up to five different colors for multiple overlapping
        % plots (this color table is best suited if the vector fields
        % are to be overlaid on top of an (FSM) image)
        colorTable={'yellow','red','green','magenta','blue'};
    end
    
    % Only up two five plots can be colored
    if nPlots>10
        nPlots=10;
    end
    
    % Change colors
    counter=nPlots/2+1;
    for i=1:2:nPlots-1
        counter=counter-1;
        set(plotHandles(i),'Color',char(colorTable(counter)));
        set(plotHandles(i+1),'Color',char(colorTable(counter)));
    end
    
elseif size(F,2)==3
    
    % Use manual scaling instead of quiver scaling which is field-dependent
    if scale~=1
        F(:,3)=scale*F(:,3);
    end

    Z=(reshape(F(:,3),length(unique(F(:,2))),length(unique(F(:,1)))))';
    h=figure;
    surf(Z);
    axis ij
    xlabel('x');
    ylabel('y');
    zlabel('z');
    view(2)
    colorbar
    
else
    
    error('Field is not a valid VECTOR or DIVERGENCE field.');
    
end