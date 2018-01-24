function cord = makiCoord2Cord(initCoord,dataProperties,forCutoff)
%MAKICOORD2CORD converts initCoord to a cord-array for the detector
%
% SYNOPSIS: cord = makiCoord2Cord(initCoord)
%
% INPUT initCoord: initCoord-cell from makiInitCoord
%       forCutoff: only use the coordinates for establishing cutoff
%
% OUTPUT cord: cord-structure used for detectSpots-subfunctions
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 02-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(dataProperties.crop)
    dataProperties.crop = zeros(2,3);
end

crop = dataProperties.crop(:,1:3);
isCrop = any(crop,1);
crop(1,~isCrop) = 1;
crop(2,~isCrop) = dataProperties.movieSize(find(~isCrop));


nTimepoints = length(initCoord);
if nargin < 3 || isempty(forCutoff)
    forCutoff = false;
end
if ~forCutoff
    fname = 'allCoordPix';
else
    fname = 'data4MMF';
end

% write sp, COM. Should anything more be necessary, check spotfind.m for
% how the information is filled in
cord(1:nTimepoints) = ...
    struct('sp',[],'COM',[]);
for t=1:nTimepoints
    
    initCoord(t).allCoordPix = initCoord(t).allCoordPix(:,1:3)  + ...
        repmat(crop(1,:)-1,initCoord(t).nSpots,1);
    
    initCoord(t).data4MMF = initCoord(t).data4MMF(:,1:3)  + ...
        repmat(crop(1,:)-1,4,1);
    
    if ~isempty(initCoord(t).nSpots)
        for i=1:size(initCoord(t).(fname),1)
%             cord(t).sp(i).cord = ...
%                 initCoord(t).(fname)(i,[2,1,3]);
 cord(t).sp(i).cord = ...
                initCoord(t).(fname)(i,:);
        end
%         cord(t).COM = ...
%             mean(initCoord(t).(fname)(:,[2,1,3]),1);
        cord(t).COM = ...
            mean(initCoord(t).(fname),1);
    end
end