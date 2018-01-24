function data = frame2FrameTrackDisplacement_oldV( tracks,single )
% EHarry October 2011
if nargin < 2 || isempty(single)
    single = 1;
end

if isstruct(tracks)
    tracksInfo = convStruct2MatNoMS(tracks);
else
    tracksInfo = tracks;
end

tracksXYZ = tracksInfo(:,1:8:end);
tracksXYZ(:,:,2) = tracksInfo(:,2:8:end);
tracksXYZ(:,:,3) = tracksInfo(:,3:8:end);
%tracksXYZ(:,:,3) = deal(0);

%tracksXYZ = tracksXYZ(:,1:10:end,:);


tracksDisp = sqrt(sum(((tracksXYZ(:,2:end,:)-tracksXYZ(:,1:end-1,:)).^2),3));

if single
    data = tracksDisp(:);
else
    data = tracksDisp;
end

end

