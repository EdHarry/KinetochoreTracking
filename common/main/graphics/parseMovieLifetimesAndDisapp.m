function[]=parseMovieLifetimesAndDisapp(movie,trackinfo,LFTmin,preTime);
% parseMovieLifetimesAndDisapp parses tracks according to their
% disappearance and superimposes the detected spots on the tiff files
% INPUT:    movie =     X x Y x N matrix containing n tiff files of size XxY
%           trackinfo = trackInfo matrix from Khuloud's tracker
%           LFTmin  =   minimum required lifetime of objects
%           preTime =   the time period before disappearance that the
%                       objects are displayed 
% 
% OUTPUT: none
% NOTE: if you want to get a movie, then run the code
% MakeQTMovie start 'nameOfMovie'
% before running this function, and afterwards run
% MakeQTMovie finish
% AND uncomment the line 47 below

% initialize parameters
[ix,iy,numf] = size(movie);
[stx,sty] = size(trackinfo);
tsiz = round(sty/8);

%determine lifetimes of all objects in trackinfo
[lifetimesMat]=findLifetimes(trackinfo);

posLFTbin1 =  find(lifetimesMat >= LFTmin );
bin1Mat = sparse(stx,numf); bin1Mat(posLFTbin1) = 1;

%determine disappearing points of all objects in trackinfo
[disappMat]=findDisappPoint(trackinfo, preTime);

trackinfoX = trackinfo(:,1:8:sty);
trackinfoY = trackinfo(:,2:8:sty);

projXb1 = trackinfoX.*bin1Mat.*disappMat;
projYb1 = trackinfoY.*bin1Mat.*disappMat;

% display parsed lifetime pits with image
for t=1:numf-1
    imshow(movie(:,:,t),[]); hold on;
    
    plot(projXb1(:,t),projYb1(:,t),'mo','LineWidth',2);
    
    pause(0.1);
    % if the series should be written into a movie, run the ocde
    % MakeQTMovie start filename
    % before this function, and uncomment the follwong line
     MakeQTMovie addaxes
    % then run after this function:
    % MakeQTMovie finish
end

end % of function



%==========================================================================
% subfunctions
%==========================================================================


function[lifetimesMat]=findLifetimes(trackinfo);

% size
[stx,sty] = size(trackinfo);
tsiz = round(sty/8);
% lifetimes matrix has number of columns corresponding to number of frames
lifetimesMat = sparse(stx,tsiz);
for i=1:stx
    % find defined x coordinates
    xcoords = find(trackinfo(i,1:8:sty)>0);
    % find defined y coordinates
    ycoords = find(trackinfo(i,2:8:sty)>0);
    % since tracks may be gap-closed, use max-min+1 instead of length
    lifetime = max ( length(xcoords), max(xcoords)-min(xcoords)+1 );
    lifetimesMat(i,xcoords) = lifetime;
end

end % of function


%==========================================================================


function[disappMat]=findDisappPoint(trackinfo, shift);

% size
[stx,sty] = size(trackinfo);
tsiz = round(sty/8);
% lifetimes matrix has number of columns corresponding to number of frames
disappMat = sparse(stx,tsiz);
for i=1:stx
    % find defined x coordinates
    xcoords = find(trackinfo(i,1:8:sty)>0);
    % find defined y coordinates
    ycoords = find(trackinfo(i,2:8:sty)>0);
    % since tracks may be gap-closed, use max-min+1 instead of length
    disappearancePoint = max ( max(xcoords), max(ycoords) );
    % ONLY COUNT AS REAL DISAPPEARANCE IF THE VLUE ISN'T THE SAME AS THE
    % LENGTH OF THE MOVIE
    if disappearancePoint<tsiz
        startPoint = max((disappearancePoint-shift),1);
        disappMat(i,startPoint:disappearancePoint) = 1;
    end
end

end % of function

