function selectIdxL = selectRois(ah)
%SELECTRECTANGLES allows the selection of multiple regions on a plot
%
% SYNOPSIS selectIdxCell = selectRois(ah)
%
% INPUT ah : (opt) axes handle. Default: gca
%
% OUTPUT selectIdxL : nPoints-by-mRois logical array indicating
%                     whether a point is inside a rectangle
%        
% c: jonas 5/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
if nargin == 0 || isempty(ah)
    ah = gca;
end

% read point coordinates
ph = findall(ah,'Type','line');
if length(ph) == 1
    % good
    X = get(ph,'XData');
    Y = get(ph,'YData');
    nData = length(X);
else
    % theoretically, it would be possible to return indices for multiple
    % plots in cell arrays, but it'll be too long to implement tonight
    error('axes can only have one plot')
end

% preassign output
selectIdxL = false(nData,100);
deltaNRoi = 100;

% remember hold state
holdOn = ishold;
hold on

% loop till the user does not want to select any more ROIs. 
nRoi = 1;
done = false;

while ~done
    
    % get polygon via getline
    [xRoi,yRoi] = getline(ah,'closed');
    
    % draw polygon
    lh=line(xRoi,yRoi,'Color',extendedColors(nRoi));
    
    % find indices of points inside polygon
    selectIdxL(:,nRoi) = inpolygon(X,Y,xRoi,yRoi);
    
    % ask for more
    cont = myQuestdlg('Do you want to select another ROI?','',...
        'Yes','No','Remove last','Remove last and quit','Yes');
    switch cont
        case 'Yes'
            % up nRoi, check whether we have enough space
            if nRoi == size(selectIdxL,2)
                % add some more columns
                selectIdxL = [selectIdxL,false(nData,deltaNRoi)];
            end
            nRoi = nRoi + 1;
            % not done
            
        case 'Remove last'
            % remove the line, redo current ROI
            delete(lh)
            
            % not done
            
        case 'Remove last and quit'
            % remove the line, exit
            delete(lh)
            selectidxL(:,nRoi:end) = [];
            done = true;
            
            
        otherwise
            % done and exit
            
            % remove superfluous cols
            selectIdxL(:,nRoi+1:end) = [];
            
            % done
            done = true;           
        
    end
end

if ~holdOn
    hold off
end

