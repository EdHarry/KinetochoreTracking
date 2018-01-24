function [pixList] = orderConnectedPixels(bwMask)
% ORDERCONNECTEDPIXELS sorts pixels that define one or more edges
%
% INPUT: bwMask: binary mask containing n curves or edges
%
% OUTPUT: pixList: n x 1 cell array where n is the number of disconnected
%                  edges found in the image. the ith entry contains
%                  segment i's pixels in order by connectivity. if two
%                  edges cross each other, the intersection is removed and
%                  four segments are returned. if an edge is closed (like a
%                  cell outline), the first pixel in the list will be the
%                  closest to the upper left corner of the image.
%
% 2008 - Kathryn


% initialize in case no segments found
segCount=1;
pixList=[];

% get pix indices of BW image
[rows,cols]=find(bwMask==1);
pixInd=find(bwMask);

while ~isempty(pixInd)

    % neighbors are 1 or sqrt(2) pixels from each other if they are 4- or
    % 8-connected, respectively
    D=createDistanceMatrix([rows cols],[rows cols]);
    D(D>sqrt(2))=0;

    % take out those with too many neighbors, find D again
    tooManyNeighborsIdx=find(sum(D>0,2)>2); % if 3 or more
    bwMask(pixInd(tooManyNeighborsIdx))=0;

    % find possible neighbors again
    [rows,cols]=find(bwMask);
    pixInd=find(bwMask);
    D=createDistanceMatrix([rows cols],[rows cols]);
    D(D>sqrt(2))=0;

    % possible endpoints for lines (only have 1 neighbor)
    segmentEnds=find(sum(D,2)<=sqrt(2));

    if isempty(segmentEnds) && ~isempty(pixInd) % then segment is closed
        segmentEnds=1;
    end

    % iterate through segment ends and find list of adjoining pixels
    while ~isempty(segmentEnds)
        % here's the next segment end to work on
        seg=segmentEnds(1);

        nextIdx=1;
        while ~isempty(nextIdx)
            cands=find(D(seg(end),:)); % indices of possible next pixels
            nextIdx=setdiff(cands,seg); % get rid of ones already in the segment
            if ~isempty(nextIdx) % there's still at least one possibility
                seg=[seg nextIdx(1)]; % so we add it to the list
            else
                % we are done with the current segment, so now we
                % remove seg from the list of segment ends
                segmentEnds=setdiff(segmentEnds,seg); 
            end

        end

        % assign pix indices to struct of segments
        pixList{segCount,1}=pixInd(seg);

        segCount=segCount+1;
    end
    
    % remove pix just recorded from BW image, then get indices of remaining
    bwMask(cat(1,pixList{:,1}))=0;
    [rows,cols]=find(bwMask==1);
    pixInd=find(bwMask);

end




