function cost = fitF2FDMultiple( dataStruct,type,aligned )
% EHarry October 2011
% data = [];
%
% for i = 1:length(tracksM)
%     dataTemp = frame2FrameTrackDisplacement(tracksM(i).tracks);
%     data(end+1:end+length(dataTemp)) = dataTemp;
% end
%
% data = data';

if nargin < 2
    type = [];
end

if nargin < 3 || isempty(aligned)
    aligned = 0;
end

data = frame2FrameTrackDisplacments_all( dataStruct,type,aligned,0 );

eval(['data = data.data_',type,';']);

[~, ~, inlierIdx, outlierIdx] = robustMean(data,[],[],1);

bic = fitFrameToFrameDisplacements(data(inlierIdx));

cost = bic + length(outlierIdx); % the cost is the bic plus the number of outliers


end

