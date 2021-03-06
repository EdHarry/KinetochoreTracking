function initCoord = makiSlist2Coord(slist,dataProperties)
% EHarry Dec 2011


if isempty(dataProperties.crop)
    dataProperties.crop = zeros(2,3);
end

crop = dataProperties.crop(:,1:3);
isCrop = any(crop,1);
crop(1,~isCrop) = 1;
crop(2,~isCrop) = dataProperties.movieSize(find(~isCrop));


pixelSize = [dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z];

nTimepoints = length(slist);
% if nargin < 3 || isempty(forCutoff)
%     forCutoff = false;
% end
% if ~forCutoff
%     fname = 'allCoordPix';
% else
%     fname = 'data4MMF';
% end

% write sp, COM. Should anything more be necessary, check spotfind.m for
% how the information is filled in
% cord(1:nTimepoints) = ...
%     struct('sp',[],'COM',[]);
initCoord(1:nTimepoints) = struct('allCoord',[],'allCoordPix',[],'nSpots',[],'amp',[]);
for t=1:nTimepoints
    
    %     initCoord(t).allCoordPix = initCoord(t).allCoordPix(:,1:3)  + ...
    %         repmat(crop(1,:)-1,initCoord(t).nSpots,1);
    %
    %     initCoord(t).data4MMF = initCoord(t).data4MMF(:,1:3)  + ...
    %         repmat(crop(1,:)-1,4,1);
    
    sp = slist(t).sp;
    
    Q = diag(slist(t).statistics.Q);
    qAmp = diag(slist(t).statistics.qAmp);
    
    initCoord(t).nSpots = length(sp);
    
    % if ~isempty(sp)
    for i=1:length(sp)
        %             cord(t).sp(i).cord = ...
        %                 initCoord(t).(fname)(i,[2,1,3]);
        initCoord(t).allCoordPix(i,1:3) = sp(i).cord  - crop(1,:);
        initCoord(t).allCoordPix(i,4:6) = Q(((i-1)*3)+1:((i-1)*3)+3)';
        initCoord(t).amp(i,1) = sp(i).amp;
        initCoord(t).amp(i,2) = qAmp(i);
    end
    
    
    initCoord(t).allCoord = initCoord(t).allCoordPix.*repmat(pixelSize,initCoord(t).nSpots,2);
    
    %end
end