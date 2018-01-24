function lineMapOut  = drawpoly(polyIn,matSize)

%%Description%%
%
% lineMapOut  = drawpoly(polyIn,matSize)
%
%Returns a matrix containing nonzeros in points which are traversed by the edges
%of the input polygon polyIn. The values in these points correlates to the
%approximate index in the input polygon which they came from.
%
%Use this output along with an image of size matSize to perform a quick and
%dirty linescan/polygonscan! 
%However, if you do this be careful not to use a
%polygon which intersects with itself as this will cause gaps in the
%resulting linescan.
%
%%Input%%
%
% polyIn >>>> This is a 2xM vector (where M > 1 ) of the polygon's vertices.
%                    Can also be a line.
%
% matSize>>>Optional. The size of the binary matrix the polygon will be
%                    returned in.
%
%
%Hunter Elliott, 3/2009

%% ------ Input ------ %%

if nargin < 1 || isempty(polyIn)
    error('Must input a polygong!')
end

[check,nPoints] = size(polyIn);


%TEMP TEMP
tmp = polyIn;
polyIn(1,:) = tmp(2,:);
polyIn(2,:) = tmp(1,:);

if check ~= 2 || nPoints < 2
    error('The input polygon must be a 2xM vector where M >= 2')
end

%Check that the polygon has only positive values for it's vertex coord.
if min(polyIn(:)) < 0
    error('Polygon vertices must have positive coordinates!')
end

if nargin < 2 || isempty(matSize)    
    matSize(1) = max(polyIn(1,:))+1;
    matSize(2) = max(polyIn(2,:))+1;
end


%Init mask
lineMapOut = zeros(matSize);

%Loop through the points and draw lines connecting them

%Over-initialize output indices
indicesOut = nan(1,nPoints*2);
currInd = 1;
for iPt = 2:nPoints
    
    %Get the indices for readability
    iM = polyIn(1,iPt-1);
    fM = polyIn(1,iPt);
    iN = polyIn(2,iPt-1);
    fN = polyIn(2,iPt);
   
    
    %Fill in the polygon edges with 1s
    
    %Check for vertical or horizontal lines
    if iM-fM == 0
        iM = floor(iM);
        iM(iM<1) = 1;
        iN = floor(iN);
        iN(iN<1)=1;
        fN = floor(fN);
        fN(fN < 1) = 1;
        nPts = (fN-iN);        
        currInd = linspace(iPt-1,iPt,nPts);
        lineMapOut(iM,iN:fN) = currInd;                        
    elseif iN-fN == 0
        iM = floor(iM);
        iM(iM<1) = 1;
        iN = floor(iN);
        iN(iN<1)=1;
        fM = floor(fM);
        fM(fM < 1) = 1;
        nPts = (fM-iM);        
        currInd = linspace(iPt-1,iPt,nPts);        
        lineMapOut(iM:fM,iN) = currInd;        
    else %Fill in diagonal lines with no gaps         
        if abs(fN-iN) > abs(fM-iM)
            tmpN = linspace(iN,fN,abs(fN-iN)+1*2); %x values. Double the points to avoid gaps from rounding.
            mLine = (fM-iM)/(fN-iN);%slope
            bLine = iM-mLine*iN;  %intercept
            tmpM = min(ceil(tmpN .* mLine + bLine),matSize(1));    %y values            
        else
            tmpM = linspace(iM,fM,abs(fM-iM)+1*2); %x values
            mLine = (fN-iN)/(fM-iM);%slope
            bLine = iN-mLine*iM;  %intercept
            tmpN = min(ceil(tmpM .* mLine + bLine),matSize(2));    %y values                        
        end
        tmpM = floor(tmpM);%Round the values and remove zeros after the line has been calculated
        tmpM(tmpM < 1) = 1;
        tmpN = floor(tmpN);
        tmpN(tmpN < 1) = 1;
        
        %Get the indices corresponding to these points;
        nPts =  length(tmpN);
        currInd = linspace(iPt-1,iPt,nPts);        
        
        for k = 1:nPts;            
            lineMapOut(tmpM(k),tmpN(k)) = currInd(k);                
        end
    end     
    
end

