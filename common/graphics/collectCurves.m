function collectionFigure = collectCurves(figureHandles, recolor, showLegend)
%COLLECTCURVES will copy curves from several figures into one figure
%
% SYNOPSIS collectCurves(figureHandles, recolor)
%
% INPUT    figureHandles : vector of figure handles
%          recolor       : (opt) [0/{1}] whether or not to recolor the curves
%          showLegend    : (opt) [0/{1}] whether or not to show legend
%
% OUTPUT   collectionFigure : handle to the figure with the collected curves
%
% REMARKS  To better identify curves, it might be helpful to tag them
%           first, e.g. by using plot(x,y,'Tag','myDescription'). After
%           collection, the functions '(un)hideErrorbars' can be very
%           convenient for clarity.
%
%
%c: jonas 04/04
% 11/07 - support for bars (ML7)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------
% test input
%------------
if nargin == 0 || isempty(figureHandles) || ~all(ishandle(figureHandles))
    error('please specify valid figure handles as input for COLLECTCURVES')
end
if nargin < 2 || isempty(recolor)
    recolor = 1;
end
if nargin < 3 || isempty(showLegend)
    showLegend = 1;
end
%-------------

%-----------------
% make new figure
newFH = figure('Name','Collected Curves');
newAxH = axes;

% reshape figureHandles so that we can run the loop correctly
figureHandles = figureHandles(:)';

lineCt = 1;
colorCt = 1;
legendCt = 1;
%----------------------


%---------------------
% loop through figure handles and copy all lines. Add a tag so that you
% will know where the curves came from
for fh = figureHandles

    % to be sure: kill legends
    legH = findall(fh,'Tag','legend');
    if ~isempty(legH)
        %figHadLegend = 1;
        delete(legH);
    else
        %figHadLegend = 0;
    end

    % find lines & reshape. Order backwards, because new line handles are
    % appended to the figure children at the top, and we want to observe
    % the sequence in which the lines were plotted
    lineHandles = findall(fh,'Type','line');
    % check whether there are bars, and whether any of the lineHandles
    % correspond to the y-axis
    barHandles = findall(fh,'Type','patch');
    if ~isempty(barHandles)
        % remove the lineHandles whose yData is [0,0]
        badLines = findall(lineHandles,'YData',[0,0]);
        hIdx = ismember(lineHandles,badLines);
        lineHandles(hIdx) = [];
    end
    
    lineHandles = lineHandles(end:-1:1)';
    if ~isempty(lineHandles)
        % loop through the lines, copy to figure and update
        for lh = lineHandles

            % copy into new figure
            newH(lineCt) = copyobj(lh,newAxH); %#ok<AGROW>

            % change tag
            oldTag = get(newH(lineCt),'Tag');
            newTag = ['fig-' num2str(fh) ' ' oldTag];
            set(newH(lineCt),'Tag',newTag);

            % change color
            if recolor
                if ~isempty(findstr(oldTag,'errorBar'))
                    % set line color
                    set(newH(lineCt),'Color',extendedColors(colorCt-1));

                elseif ~isempty(findstr(oldTag, 'TAfit'))
                    % don't change color

                else
                    % set new color
                    set(newH(lineCt),'Color',extendedColors(colorCt));
                    colorCt = colorCt+1;
                end
            end

            % collect tags for legend
            if ~isempty(findstr(oldTag,'errorBar')) || ~isempty(findstr(oldTag,'TAfit'))
                % do not add to legend
            else
                legendCell{legendCt,1} = newTag; %#ok<AGROW>
                legendLineH(legendCt,1) = newH(lineCt); %#ok<AGROW>
                legendCt = legendCt + 1;
            end

            lineCt = lineCt + 1;
        end
    end
    %
    
    % checked for barHandles above
    barHandles = barHandles(end:-1:1);
    if ~isempty(barHandles)
        for bh = barHandles

            % copy into new figure
            newH(lineCt) = copyobj(bh,newAxH); %#ok<AGROW>

            % change tag
            oldTag = get(newH(lineCt),'Tag');
            newTag = ['fig-' num2str(fh) ' ' oldTag];
            set(newH(lineCt),'Tag',newTag);

            % change color
            if recolor
                if ~isempty(findstr(oldTag,'errorBar'))
                    % set line color
                    set(newH(lineCt),'Color',extendedColors(colorCt-1));

                elseif ~isempty(findstr(oldTag, 'TAfit'))
                    % don't change color

                else
                    % set new color
                    set(newH(lineCt),'FaceColor',extendedColors(colorCt));
                    set(newH(lineCt),'EdgeColor',extendedColors(colorCt));
                    colorCt = colorCt+1;
                end
            end

            % collect tags for legend
            if ~isempty(findstr(oldTag,'errorBar')) || ~isempty(findstr(oldTag,'TAfit'))
                % do not add to legend
            else
                legendCell{legendCt,1} = newTag; %#ok<AGROW>
                legendLineH(legendCt,1) = newH(lineCt); %#ok<AGROW>
                legendCt = legendCt + 1;
            end

            lineCt = lineCt + 1;
        end
    end

    %     % turn legend back on - does not work for some reason
    %     if figHadLegend
    %         axH = findall(fh,'Type','axes');
    %         legend(axH,'show');
    %     end
end

%---------------------
% now show the legend.
if showLegend
    lh=legend(legendLineH,legendCell);
    lHandles = get(lh);
    if isfield(lHandles,'Interpreter')
        set(lh,'Interpreter','none')
    end
end

%-----------------------
% nargout if asked for

if nargout > 0
    collectionFigure = newFH;
end