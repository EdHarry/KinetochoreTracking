function initCoord = convertMovieInfo2InitCoord( movieInfo,dataProperties )
%CONVERTMOVIEINFO2INITCOORD convery coords from a 'movieInfo' struct to an
% 'initCoord' struct
%   EHarry March 2012

numTimePoints = length(movieInfo);

pixelSize = [dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z];

initCoord(1:numTimePoints) = struct('allCoord',[],'allCoordPix',[],'nSpots',[],'amp',[]);

for t = 1:numTimePoints
    
    if ~isempty(movieInfo(t).xCoord)
        
        initCoord(t).allCoordPix = [movieInfo(t).xCoord(:,1) movieInfo(t).yCoord(:,1) movieInfo(t).zCoord(:,1) movieInfo(t).xCoord(:,2) movieInfo(t).yCoord(:,2) movieInfo(t).zCoord(:,2)];
        
        initCoord(t).nSpots = size(movieInfo(t).xCoord,1);
        
        initCoord(t).allCoord = initCoord(t).allCoordPix .* repmat(pixelSize,initCoord(t).nSpots,2);
        
        initCoord(t).amp = movieInfo(t).amp;
        
    end
    
end
end

