function dataStruct = KTPairTracking(dataStruct)
% improves the tracks of KT pairs
%   EHarry

%% initialise
% end if no updatedClass or sisterList or frameAlignment then stop
if isempty(dataStruct.sisterList) || isempty(dataStruct.sisterList(1).trackPairs) || isempty(dataStruct.updatedClass) || isempty(dataStruct.frameAlignment)
    return
end

% info for loading movie
movieFileName = fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName);
% use deconvolved image or not
if isfield(dataStruct.dataProperties,'decon') && ~isempty(dataStruct.dataProperties.decon)
    decon = dataStruct.dataProperties.decon;
else
    decon = 0;
end
% get cropping if any
if isfield(dataStruct.dataProperties,'crop')
    crop = dataStruct.dataProperties.crop;
else
    crop = [];
end
% get psf
if isfield(dataStruct.dataProperties,'psfSigma') && ~isempty(dataStruct.dataProperties.psfSigma)
    psfSigma = dataStruct.dataProperties.psfSigma;
else
    %get initial guess of PSF sigma
    filterPrm = dataStruct.dataProperties.FILTERPRM;
    psfSigma = filterPrm([1 3]);
    clear filterPrm
end
%get camera bit depth
if ~isfield(dataStruct.dataProperties,'bitDepth') || isempty(dataStruct.dataProperties.bitDepth)
    bitDepth = 16;
else
    bitDepth = dataStruct.dataProperties.bitDepth;
end
% get pixelSize
pixelSize = [dataStruct.dataProperties.PIXELSIZE_XY dataStruct.dataProperties.PIXELSIZE_XY dataStruct.dataProperties.PIXELSIZE_Z];
% get full psf
psfSigma = psfSigma([1 1 2]);

% get frame alignment for rotated coords and updatedClss for unaligned
% indexes and initCoord for re-fitting
frameAlignment = dataStruct.frameAlignment;
updatedClass = dataStruct.updatedClass;
initCoord = dataStruct.initCoord;
planeFit = dataStruct.planeFit;

% get no of time points
nFrames = dataStruct.dataProperties.movieSize(4);

% load
movie = readOMEMatFile(movieFileName,1:nFrames,1,decon,crop);
movie = movie ./ (2^bitDepth - 1);
% get pixel size
[xSize,ySize,zSize,~] = size(movie);

% get max sister distance
maxSisDis = dataStruct.dataProperties.groupSisters.maxDist;

% sudo work, transform max average sister dis to max possible sister dis
maxSisDis = 4.375*maxSisDis^2 - 10.62*maxSisDis + 7.894;

% get kinetochore search radii
centreSR = dataStruct.dataProperties.tracksParam.costMatrices(1).parameters.maxSearchRadius(1:2) / 2;

% get cell phase
phase = catStruct(1,'dataStruct.updatedClass.phase',[],1); % phases of the cell cycle
firstFrameAPhase = find(phase=='a',1); % find the first frame theat the cell goes into aphase, if it does at all
if isempty(firstFrameAPhase) % if the cell did not go into aphase
    endMetaphase = nFrames; % test tracks until the end of the movie
else
    endMetaphase = firstFrameAPhase - 1; % else test upto the time just before aphase
end

% get tracks
tracks = dataStruct.tracks;

% get sister track indicies
sList = dataStruct.sisterList(1).trackPairs(:,1:2);
nSis = size(sList,1);

%% improve anaphase detection
% improve anaphase detection by looking for deviations in average sister
% separation

% onlt if there was a planeFit

if ~isempty(planeFit)
    
    % generate align sisterList
    sisterList = makiConstructAlignedSisters_allCoords(dataStruct);
    [sisterList.coords1] = sisterList.coords1Aligned;
    [sisterList.coords2] = sisterList.coords2Aligned;
    
    avSep = zeros(endMetaphase,1);
    for t = 1:endMetaphase
        tmp = NaN(nSis,1);
        for iSis = 1:nSis
            coords1 = sisterList(iSis).coords1(t,1:3);
            coords2 = sisterList(iSis).coords2(t,1:3);
            tmp(iSis) = normList(coords1-coords2);
        end
        avSep(t) = robustMean(tmp);
    end
    
    % take mean
    avSepM = robustMean(avSep);
    
    % do robust exp fitting
    [~,~,idx] = robustExponentialFit2(avSep);
    idx = setxor(1:endMetaphase,idx);
    
    % go though outliers, if they are at the end of the movie and their sep is
    % above average then make them anaphase frames
    newA = [];
    if ~isempty(idx)
        loop = 1;
        count = -1;
        while loop && count < endMetaphase
            count = count + 1;
            if idx(end) == endMetaphase - count && avSep(idx(end)) > avSepM
                newA = [idx(end); newA];%#ok<AGROW>
                idx(end) = [];
            else
                loop = 0;
            end
        end
    end
    
else
    newA = [];
end

% get times before and after A phase
timesBeforeAphase = setxor(1:endMetaphase,newA)';
timesAfterAPhase = setxor(1:nFrames,timesBeforeAphase)';
if ~isempty(timesAfterAPhase)
    timesAfterAPhase = returnRightVector(timesAfterAPhase,[],'r');
else
    timesAfterAPhase = [];
end

% initialise index matrix
indexMat = [];


%% cluster sisters
% get coords and detect outliers
for iSis = 1:nSis
    % get coords and amps
    coords1 = getCoords(tracks(sList(iSis,1)),nFrames);
    coords1(:,4) = getAmp(tracks(sList(iSis,1)),nFrames);
    coords2 = getCoords(tracks(sList(iSis,2)),nFrames);
    coords2(:,4) = getAmp(tracks(sList(iSis,2)),nFrames);
    
    % detect outliers
    outliers1 = 0;
    outliers2 = 0;
    while ~isempty(outliers1) || ~isempty(outliers2)
        [outliers1,outliers2] = detectOuliers(coords1,coords2,timesBeforeAphase,maxSisDis);
        
        % remove outliers (make NaN)
        coords1(outliers1,:) = NaN;
        coords2(outliers2,:) = NaN;
    end
    
    % find not NaNs on either sister
    coordIdx = ~isnan(coords1(:,1)) | ~isnan(coords2(:,1));
    
    % cluster coords
    clusters = bwconncomp([coordIdx coordIdx]);
    
    % get original indexes
    featIdx1 = getFeatIdx(tracks(sList(iSis,1)),nFrames)';
    featIdx2 = getFeatIdx(tracks(sList(iSis,2)),nFrames)';
    
    % cluster indicies and add to matrix
    for iCluster = 1:clusters.NumObjects
        idx = clusters.PixelIdxList{iCluster}(1:end/2);
        idx = setxor(1:nFrames,idx);
        featIdx1_tmp = featIdx1;
        featIdx1_tmp(idx) = NaN;
        featIdx2_tmp = featIdx2;
        featIdx2_tmp(idx) = NaN;
        indexMat = [indexMat; featIdx1_tmp; featIdx2_tmp];%#ok<AGROW>
    end
    
end

% get no. clustered sisters
nSis = size(indexMat,1)/2;

% this will be a list of previous extra found sisters incase the program
% gets stuck in a loop
previousExtraSis = struct('pExtraSis',[]);

% this will be a list of previous potential new spot locations to avoid
% repeats
previousNewSpots(1:nFrames) = struct('coords',[]);

% start fittingloop 2
fittingLoop2 = true;

while fittingLoop2
    
    % start fitting loop 1
    fittingLoop1 = true;
    
    while fittingLoop1
        %% STEP 1: FILL IN GAPS IN SINGLE SISTER TRACKS
        
        disp('')
        progressText(0);
        
        % forward in time
        for iFrame = 1:nFrames-1
            % assign target frame
            targetFrame = iFrame + 1;
            
            % sister to attempt to fit (sisIdx,{1,2})
            sisToFit = [];
            for iSis = 1:nSis
                idx = (iSis-1)*2 + 1;
                if (~isnan(indexMat(idx,iFrame)) && ~isnan(indexMat(idx,targetFrame)) && ~isnan(indexMat(idx+1,iFrame)) && isnan(indexMat(idx+1,targetFrame)))
                    sisToFit = [sisToFit; iSis 2];%#ok<AGROW>
                elseif (~isnan(indexMat(idx,iFrame)) && isnan(indexMat(idx,targetFrame)) && ~isnan(indexMat(idx+1,iFrame)) && ~isnan(indexMat(idx+1,targetFrame)))
                    sisToFit = [sisToFit; iSis 1];%#ok<AGROW>
                end
            end
            
            % if no sisters to fit in this frame (ok, next frame actually then continue to next frame)
            if isempty(sisToFit)
                progressText(iFrame/(nFrames-1),'forward processing single missing timepoints')
                continue
            end
            
            % first collect indicies of coords not to condider
            toRem = zeros(size(sisToFit,1),1);
            for idx = 1:size(sisToFit,1)
                iSis = sisToFit(idx,1);
                if sisToFit(idx,2) == 1
                    toRem(idx) = indexMat((iSis-1)*2 + 2,targetFrame);
                else
                    toRem(idx) = indexMat((iSis-1)*2 + 1,targetFrame);
                end
            end
            
            for iSis = 1:nSis
                idx = (iSis-1)*2 + 1;
                if (~isnan(indexMat(idx,iFrame)) && ~isnan(indexMat(idx+1,iFrame)) && ~isnan(indexMat(idx,targetFrame)) && ~isnan(indexMat(idx+1,targetFrame)))
                    toRem = [toRem; indexMat(idx,targetFrame); indexMat(idx+1,targetFrame)];%#ok<AGROW>
                end
            end
            
            % get coords in the target frame and this frame
            coords = frameAlignment(targetFrame).alignedCoord(:,1:3);
            coordsThisFrame = frameAlignment(iFrame).alignedCoord(:,1:3);
            % get index to coords to consider for fitting
            toConsider = setxor(1:size(coords,1),toRem);
            
            % crop coords
            coords2 = coords(toConsider,:);
            
            % make coords1 which is the expected positions of sisters being tracked
            % and current dis which is the current sister dis in this frame
            % also store the current angle with the normal
            coords1 = zeros(size(sisToFit,1),3);
            currentDis = zeros(size(sisToFit,1),1);
            currentAngle = zeros(size(sisToFit,1),1);
            for idx = 1:size(sisToFit,1)
                iSis = sisToFit(idx,1);
                % make displacement vector
                if sisToFit(idx,2) == 1
                    colIdx = (iSis-1)*2 + 2;
                    otherIdx = (iSis-1)*2 + 1;
                    disVec = coords(indexMat(colIdx,targetFrame),:) - coordsThisFrame(indexMat(colIdx,iFrame),:);
                else
                    colIdx = (iSis-1)*2 + 1;
                    otherIdx = (iSis-1)*2 + 2;
                    disVec = coords(indexMat(colIdx,targetFrame),:) - coordsThisFrame(indexMat(colIdx,iFrame),:);
                end
                % reverse x direction if this is an anaphase frame
                if ismember(targetFrame,timesAfterAPhase)
                    disVec(1) = -disVec(1);
                end
                % make expected coords
                coords1(idx,:) = coordsThisFrame(indexMat(otherIdx,iFrame),:) + disVec;
                % get current sister dis
                dis = normList(coordsThisFrame(indexMat(otherIdx,iFrame),:) - coordsThisFrame(indexMat(colIdx,iFrame),:));
                vecN = (coordsThisFrame(indexMat(otherIdx,iFrame),1) - coordsThisFrame(indexMat(colIdx,iFrame),1)) ./ dis;
                currentDis(idx) = dis;
                currentAngle(idx) = acos(abs(vecN(1)));
            end
            
            % make lap cost mat
            costMat = -ones(size(coords1,1),size(coords2,1));
            % loop paris
            for i = 1:size(coords1,1)
                for j = 1:size(coords2,1)
                    % get dis between sisters in traget frame
                    iSis = sisToFit(i,1);
                    if sisToFit(i,2) == 1
                        colIdx = (iSis-1)*2 + 2;
                        newDis = normList(coords(indexMat(colIdx,targetFrame),:)-coords2(j,:));
                        vecN = (coords(indexMat(colIdx,targetFrame),1)-coords2(j,1)) ./ newDis;
                        newAngle = acos(abs(vecN(1)));
                    else
                        colIdx = (iSis-1)*2 + 1;
                        newDis = normList(coords(indexMat(colIdx,targetFrame),:)-coords2(j,:));
                        vecN = (coords(indexMat(colIdx,targetFrame),1)-coords2(j,1)) ./ newDis;
                        newAngle = acos(abs(vecN(1)));
                    end
                    dis = normList(coords1(i,:)-coords2(j,:));
                    % if a) dis between sister is no too large and if the distance
                    % from the expected position isn't too large then consider,
                    % also i the cnage in angle with normal is not too large
                    % only have maxSisDis constraint in non-anaphase frames
                    if ismember(targetFrame,timesAfterAPhase)
                        sisDisConstraint = Inf;
                    else
                        sisDisConstraint = maxSisDis;
                    end
                    if  newDis < sisDisConstraint && dis < 0.5 && abs(newDis-currentDis(i)) < 0.4 && abs(newAngle - currentAngle(i)) < 20*pi/180
                        % generate cost
                        costMat(i,j) = dis * abs(1 - (newDis / currentDis(i)));
                    end
                end
            end
            
            % if no pairs possible then continue to next frame
            if all(costMat(:) == -1)
                progressText(iFrame/(nFrames-1),'forward processing single missing timepoints')
                continue
            end
            
            % perform lap and corrent links
            links = lap(costMat,[],[],1);
            links = double(links(1:size(coords1,1)));
            links(links > size(coords2,1)) = NaN;
            
            %     if any(~isnan(links))
            %         disp('new links!');
            %     end
            
            % update links to coords that were considered
            links(~isnan(links)) = toConsider(links(~isnan(links)));
            
            % add to index matrix
            for idx = 1:size(sisToFit,1)
                iSis = sisToFit(idx,1);
                sis = sisToFit(idx,2);
                indexMat((iSis-1)*2 + sis,targetFrame) = links(idx);%#ok<AGROW>
            end
            progressText(iFrame/(nFrames-1),'forward processing single missing timepoints')
        end
        
        % backwards in time
        disp('')
        progressText(0);
        
        for iFrame = nFrames:-1:2
            % assign target frame
            targetFrame = iFrame - 1;
            
            % sister to attempt to fit (sisIdx,{1,2})
            sisToFit = [];
            for iSis = 1:nSis
                idx = (iSis-1)*2 + 1;
                if (~isnan(indexMat(idx,iFrame)) && ~isnan(indexMat(idx,targetFrame)) && ~isnan(indexMat(idx+1,iFrame)) && isnan(indexMat(idx+1,targetFrame)))
                    sisToFit = [sisToFit; iSis 2];%#ok<AGROW>
                elseif (~isnan(indexMat(idx,iFrame)) && isnan(indexMat(idx,targetFrame)) && ~isnan(indexMat(idx+1,iFrame)) && ~isnan(indexMat(idx+1,targetFrame)))
                    sisToFit = [sisToFit; iSis 1];%#ok<AGROW>
                end
            end
            
            % if no sisters to fit in this frame (ok, next frame actually then continue to next frame)
            if isempty(sisToFit)
                progressText((nFrames-iFrame+1)/(nFrames-1),'backward processing single missing timepoints')
                continue
            end
            
            % first collect indicies of coords not to condider
            toRem = zeros(size(sisToFit,1),1);
            for idx = 1:size(sisToFit,1)
                iSis = sisToFit(idx,1);
                if sisToFit(idx,2) == 1
                    toRem(idx) = indexMat((iSis-1)*2 + 2,targetFrame);
                else
                    toRem(idx) = indexMat((iSis-1)*2 + 1,targetFrame);
                end
            end
            
            for iSis = 1:nSis
                idx = (iSis-1)*2 + 1;
                if (~isnan(indexMat(idx,iFrame)) && ~isnan(indexMat(idx+1,iFrame)) && ~isnan(indexMat(idx,targetFrame)) && ~isnan(indexMat(idx+1,targetFrame)))
                    toRem = [toRem; indexMat(idx,targetFrame); indexMat(idx+1,targetFrame)];%#ok<AGROW>
                end
            end
            
            % get coords in the target frame and this frame
            coords = frameAlignment(targetFrame).alignedCoord(:,1:3);
            coordsThisFrame = frameAlignment(iFrame).alignedCoord(:,1:3);
            % get index to coords to consider for fitting
            toConsider = setxor(1:size(coords,1),toRem);
            
            % crop coords
            coords2 = coords(toConsider,:);
            
            % make coords1 which is the expected positions of sisters being tracked
            % and current dis which is the current sister dis in this frame
            % also store the current angle with the normal
            coords1 = zeros(size(sisToFit,1),3);
            currentDis = zeros(size(sisToFit,1),1);
            currentAngle = zeros(size(sisToFit,1),1);
            for idx = 1:size(sisToFit,1)
                iSis = sisToFit(idx,1);
                % make displacement vector
                if sisToFit(idx,2) == 1
                    colIdx = (iSis-1)*2 + 2;
                    otherIdx = (iSis-1)*2 + 1;
                    disVec = coords(indexMat(colIdx,targetFrame),:) - coordsThisFrame(indexMat(colIdx,iFrame),:);
                else
                    colIdx = (iSis-1)*2 + 1;
                    otherIdx = (iSis-1)*2 + 2;
                    disVec = coords(indexMat(colIdx,targetFrame),:) - coordsThisFrame(indexMat(colIdx,iFrame),:);
                end
                % reverse x direction if this is an anaphase frame
                if ismember(targetFrame,timesAfterAPhase)
                    disVec(1) = -disVec(1);
                end
                % make expected coords
                coords1(idx,:) = coordsThisFrame(indexMat(otherIdx,iFrame),:) + disVec;
                % get current sister dis
                dis = normList(coordsThisFrame(indexMat(otherIdx,iFrame),:) - coordsThisFrame(indexMat(colIdx,iFrame),:));
                vecN = (coordsThisFrame(indexMat(otherIdx,iFrame),1) - coordsThisFrame(indexMat(colIdx,iFrame),1)) ./ dis;
                currentDis(idx) = dis;
                currentAngle(idx) = acos(abs(vecN(1)));
            end
            
            % make lap cost mat
            costMat = -ones(size(coords1,1),size(coords2,1));
            % loop paris
            for i = 1:size(coords1,1)
                for j = 1:size(coords2,1)
                    % get dis between sisters in traget frame
                    iSis = sisToFit(i,1);
                    if sisToFit(i,2) == 1
                        colIdx = (iSis-1)*2 + 2;
                        newDis = normList(coords(indexMat(colIdx,targetFrame),:)-coords2(j,:));
                        vecN = (coords(indexMat(colIdx,targetFrame),1)-coords2(j,1)) ./ newDis;
                        newAngle = acos(abs(vecN(1)));
                    else
                        colIdx = (iSis-1)*2 + 1;
                        newDis = normList(coords(indexMat(colIdx,targetFrame),:)-coords2(j,:));
                        vecN = (coords(indexMat(colIdx,targetFrame),1)-coords2(j,1)) ./ newDis;
                        newAngle = acos(abs(vecN(1)));
                    end
                    dis = normList(coords1(i,:)-coords2(j,:));
                    % if a) dis between sister is no too large and if the distance
                    % from the expected position isn't too large then consider,
                    % also i the cnage in angle with normal is not too large
                    if ismember(targetFrame,timesAfterAPhase)
                        sisDisConstraint = Inf;
                    else
                        sisDisConstraint = maxSisDis;
                    end
                    if  newDis < sisDisConstraint && dis < 0.5 && abs(newDis-currentDis(i)) < 0.4 && abs(newAngle - currentAngle(i)) < 20*pi/180
                        % generate cost
                        costMat(i,j) = dis * abs(1 - (newDis / currentDis(i)));
                    end
                end
            end
            
            % if no pairs possible then continue to next frame
            if all(costMat(:) == -1)
                progressText((iFrame-1)/(nFrames-1),'backward processing single missing timepoints')
                continue
            end
            
            % perform lap and corrent links
            links = lap(costMat,[],[],1);
            links = double(links(1:size(coords1,1)));
            links(links > size(coords2,1)) = NaN;
            
            %     if any(~isnan(links))
            %         disp('new links!');
            %     end
            
            % update links to coords that were considered
            links(~isnan(links)) = toConsider(links(~isnan(links)));
            
            % add to index matrix
            for idx = 1:size(sisToFit,1)
                iSis = sisToFit(idx,1);
                sis = sisToFit(idx,2);
                indexMat((iSis-1)*2 + sis,targetFrame) = links(idx);
            end
            progressText((iFrame-1)/(nFrames-1),'backward processing single missing timepoints')
        end
        
        
        %% STEP 2: FILL IN DOUBLE SISTER GAPS
        
        % forward in time
        disp('')
        progressText(0);
        for iFrame = 1:nFrames-1
            % assign target frame
            targetFrame = iFrame + 1;
            
            % sister to attempt to fit
            sisToFit = [];
            for iSis = 1:nSis
                idx = (iSis-1)*2 + 1;
                if (~isnan(indexMat(idx,iFrame)) && ~isnan(indexMat(idx+1,iFrame)) && isnan(indexMat(idx,targetFrame)) && isnan(indexMat(idx+1,targetFrame)))
                    sisToFit = [sisToFit; iSis];%#ok<AGROW>
                end
            end
            
            % if no sisters to fit in this frame (ok, next frame actually then continue to next frame)
            if isempty(sisToFit)
                progressText((iFrame)/(nFrames-1),'forward processing double missing timepoints')
                continue
            end
            
            % collect indicies of spots not to consider, tracks of sisters that
            % don't start in the traget frame
            toRem = [];
            for iSis = 1:nSis
                idx = (iSis-1)*2 + 1;
                if (~isnan(indexMat(idx,iFrame)) && ~isnan(indexMat(idx+1,iFrame)) && ~isnan(indexMat(idx,targetFrame)) && ~isnan(indexMat(idx+1,targetFrame)))
                    toRem = [toRem; indexMat(idx,targetFrame); indexMat(idx+1,targetFrame)];%#ok<AGROW>
                end
            end
            
            % get coords to consider
            if isempty(frameAlignment(targetFrame).alignedCoord)
                progressText((iFrame)/(nFrames-1),'forward processing double missing timepoints')
                continue
            else
                coordsTargetFrame = frameAlignment(targetFrame).alignedCoord(:,1:3);
            end
            if isempty(frameAlignment(iFrame).alignedCoord)
                progressText((iFrame)/(nFrames-1),'forward processing double missing timepoints')
                continue
            else
                coordsThisFrame = frameAlignment(iFrame).alignedCoord(:,1:3);
            end
            toConsider = setxor(1:size(coordsTargetFrame,1),toRem);
            coords = coordsTargetFrame(toConsider,:);
            
            % form possible pairs by cutoff in dis mat
            % onlt have max sister dis constraint in non-anaphase frames
            if ismember(targetFrame,timesAfterAPhase)
                sisDisConstraint = Inf;
            else
                sisDisConstraint = maxSisDis;
            end
            disMat = createSparseDistanceMatrix(coords,coords,sisDisConstraint,-1);
            [rows,cols] = find(disMat>0);
            pairs = [rows cols];
            pairs = sort(pairs,2);
            pairs = unique(pairs,'rows');
            
            % if no pairs possible then continue to next frame
            if isempty(pairs)
                progressText((iFrame)/(nFrames-1),'forward processing double missing timepoints')
                continue
            end
            
            % allow larger search radii for unaligned coords, first move unalinged
            % target pairs to the end of the list
            pairsUnaligned = ismember(toConsider(pairs(:,1)),updatedClass(targetFrame).unalignedIdx) | ismember(toConsider(pairs(:,2)),updatedClass(targetFrame).unalignedIdx);
            pairs = [pairs(~pairsUnaligned,:); pairs(pairsUnaligned,:)];
            nPairsAligned = sum(~pairsUnaligned);
            
            % make possible centre coords of targets, and inter sister distance and
            % angles with normals
            centreCoords2 = coords(pairs(:,1),:);
            centreCoords2 = cat(3,centreCoords2,coords(pairs(:,2),:));
            [newSisDis,vecN] = normList(centreCoords2(:,:,1) - centreCoords2(:,:,2));
            newAngles = acos(abs(vecN(:,1)));
            centreCoords2 = mean(centreCoords2,3);
            
            % make centre coords of sisters, and current inter sister dis
            centreCoords1 = zeros(size(sisToFit,1),3);
            sisDis = zeros(size(sisToFit,1),1);
            sisAngles = zeros(size(sisToFit,1),1);
            unalignedSis = false(size(sisToFit,1),1);
            for idx = 1:size(sisToFit,1)
                iSis = sisToFit(idx);
                sisIdx = indexMat((iSis-1)*2+1:(iSis-1)*2+2,iFrame);
                unalignedSis(idx) = ismember(sisIdx(1),updatedClass(iFrame).unalignedIdx) || ismember(sisIdx(2),updatedClass(iFrame).unalignedIdx);
                centreCoords1(idx,:) = mean([coordsThisFrame(sisIdx(1),:); coordsThisFrame(sisIdx(2),:)]);
                dis = normList(coordsThisFrame(sisIdx(1),:) - coordsThisFrame(sisIdx(2),:));
                vecN = (coordsThisFrame(sisIdx(1),1) - coordsThisFrame(sisIdx(2),1)) ./ dis;
                sisDis(idx) = dis;
                sisAngles(idx) = acos(abs(vecN(1)));
            end
            
            % sort unaligned pairs
            sisToFit = [sisToFit(~unalignedSis); sisToFit(unalignedSis)];
            centreCoords1 = [centreCoords1(~unalignedSis,:); centreCoords1(unalignedSis,:)];
            sisDis = [sisDis(~unalignedSis); sisDis(unalignedSis)];
            sisAngles = [sisAngles(~unalignedSis); sisAngles(unalignedSis)];
            nSisAligned = sum(~unalignedSis);
            
            % make cost mat with differnt cutoffs for aligned and unaligned pairs
            costMatLT = createSparseDistanceMatrix(centreCoords1(1:nSisAligned,:),centreCoords2(1:nPairsAligned,:),centreSR(1));
            costMatLB = createSparseDistanceMatrix(centreCoords1(nSisAligned+1:end,:),centreCoords2(1:nPairsAligned,:),centreSR(2));
            costMatRT = createSparseDistanceMatrix(centreCoords1(1:nSisAligned,:),centreCoords2(nPairsAligned+1:end,:),centreSR(2));
            costMatRB = createSparseDistanceMatrix(centreCoords1(nSisAligned+1:end,:),centreCoords2(nPairsAligned+1:end,:),centreSR(2));
            costMat = [costMatLT costMatRT; costMatLB costMatRB];
            costMat = full(costMat);
            costMat(costMat==0) = -1;
            
            % if no links possible continue to next iteration
            if all(costMat(:) == -1)
                progressText((iFrame)/(nFrames-1),'forward processing double missing timepoints')
                continue
            end
            
            % update cost
            [idxRow,idxCol] = find(costMat ~= -1);
            if size(idxRow,1) == 1
                idxRow = idxRow';
            end
            if size(idxCol,1) == 1
                idxCol = idxCol';
            end
            extraCost = abs(1 - (sisDis(idxRow) ./ newSisDis(idxCol)));
            constraint1 = abs(sisDis(idxRow) - newSisDis(idxCol)) < 0.4;
            constraint2 = abs(sisAngles(idxRow) - newAngles(idxCol)) < 20*pi/180;
            for idx = 1:length(idxRow)
                row = idxRow(idx);
                col = idxCol(idx);
                costMat(row,col) = costMat(row,col) * extraCost(idx) * constraint1(idx) * constraint2(idx);
            end
            costMat(costMat==0) = -1;
            
            % if no links possible continue to next iteration
            if all(costMat(:) == -1)
                progressText((iFrame)/(nFrames-1),'forward processing double missing timepoints')
                continue
            end
            
            % perform lap
            links = lap(costMat,[],[],1);
            links = double(links(1:size(centreCoords1,1)));
            links(links > size(centreCoords2,1)) = NaN;
            
            % get conections
            for idx = 1:size(sisToFit,1)
                iSis = sisToFit(idx);
                if ~isnan(links(idx))
                    pair = pairs(links(idx),:);
                    pair = toConsider(pair)';
                    % pair up sisters by lap
                    idxThisFrame = indexMat((iSis-1)*2+1:(iSis-1)*2+2,iFrame);
                    sisCoordsThisFrame = coordsThisFrame(idxThisFrame,:);
                    sisCoordTargetFrame = coordsTargetFrame(pair,:);
                    disMat = createDistanceMatrix(sisCoordsThisFrame,sisCoordTargetFrame);
                    pairLink = double(lap(disMat));
                    pair = pair(pairLink);
                    % fill in matrix
                    indexMat((iSis-1)*2+1:(iSis-1)*2+2,targetFrame) = pair;%#ok<AGROW>
                end
            end
            progressText((iFrame)/(nFrames-1),'forward processing double missing timepoints')
        end
        
        % backward in time
        disp('')
        progressText(0);
        for iFrame = nFrames:-1:2
            % assign target frame
            targetFrame = iFrame - 1;
            
            % sister to attempt to fit
            sisToFit = [];
            for iSis = 1:nSis
                idx = (iSis-1)*2 + 1;
                if (~isnan(indexMat(idx,iFrame)) && ~isnan(indexMat(idx+1,iFrame)) && isnan(indexMat(idx,targetFrame)) && isnan(indexMat(idx+1,targetFrame)))
                    sisToFit = [sisToFit; iSis];%#ok<AGROW>
                end
            end
            
            % if no sisters to fit in this frame (ok, next frame actually then continue to next frame)
            if isempty(sisToFit)
                progressText((nFrames-iFrame+1)/(nFrames-1),'backward processing double missing timepoints')
                continue
            end
            
            % collect indicies of spots not to consider, tracks of sisters that
            % don't start in the traget frame
            toRem = [];
            for iSis = 1:nSis
                idx = (iSis-1)*2 + 1;
                if (~isnan(indexMat(idx,iFrame)) && ~isnan(indexMat(idx+1,iFrame)) && ~isnan(indexMat(idx,targetFrame)) && ~isnan(indexMat(idx+1,targetFrame)))
                    toRem = [toRem; indexMat(idx,targetFrame); indexMat(idx+1,targetFrame)];%#ok<AGROW>
                end
            end
            
            % get coords to consider
            if isempty(frameAlignment(targetFrame).alignedCoord)
                progressText((iFrame)/(nFrames-1),'backward processing double missing timepoints')
                continue
            else
                coordsTargetFrame = frameAlignment(targetFrame).alignedCoord(:,1:3);
            end
            if isempty(frameAlignment(iFrame).alignedCoord)
                progressText((iFrame)/(nFrames-1),'backward processing double missing timepoints')
                continue
            else
                coordsThisFrame = frameAlignment(iFrame).alignedCoord(:,1:3);
            end
            toConsider = setxor(1:size(coordsTargetFrame,1),toRem);
            coords = coordsTargetFrame(toConsider,:);
            
            % form possible pairs by cutoff in dis mat
            if ismember(targetFrame,timesAfterAPhase)
                sisDisConstraint = Inf;
            else
                sisDisConstraint = maxSisDis;
            end
            disMat = createSparseDistanceMatrix(coords,coords,sisDisConstraint,-1);
            [rows,cols] = find(disMat>0);
            pairs = [rows cols];
            pairs = sort(pairs,2);
            pairs = unique(pairs,'rows');
            
            % if no pairs possible then continue to next frame
            if isempty(pairs)
                progressText((nFrames-iFrame+1)/(nFrames-1),'backward processing double missing timepoints')
                continue
            end
            
            % allow larger search radii for unaligned coords, first move unalinged
            % target pairs to the end of the list
            pairsUnaligned = ismember(toConsider(pairs(:,1)),updatedClass(targetFrame).unalignedIdx) | ismember(toConsider(pairs(:,2)),updatedClass(targetFrame).unalignedIdx);
            pairs = [pairs(~pairsUnaligned,:); pairs(pairsUnaligned,:)];
            nPairsAligned = sum(~pairsUnaligned);
            
            % make possible centre coords of targets, and inter sister distance and
            % angles with normals
            centreCoords2 = coords(pairs(:,1),:);
            centreCoords2 = cat(3,centreCoords2,coords(pairs(:,2),:));
            [newSisDis,vecN] = normList(centreCoords2(:,:,1) - centreCoords2(:,:,2));
            newAngles = acos(abs(vecN(:,1)));
            centreCoords2 = mean(centreCoords2,3);
            
            % make centre coords of sisters, and current inter sister dis
            centreCoords1 = zeros(size(sisToFit,1),3);
            sisDis = zeros(size(sisToFit,1),1);
            sisAngles = zeros(size(sisToFit,1),1);
            unalignedSis = false(size(sisToFit,1),1);
            for idx = 1:size(sisToFit,1)
                iSis = sisToFit(idx);
                sisIdx = indexMat((iSis-1)*2+1:(iSis-1)*2+2,iFrame);
                unalignedSis(idx) = ismember(sisIdx(1),updatedClass(iFrame).unalignedIdx) || ismember(sisIdx(2),updatedClass(iFrame).unalignedIdx);
                centreCoords1(idx,:) = mean([coordsThisFrame(sisIdx(1),:); coordsThisFrame(sisIdx(2),:)]);
                dis = normList(coordsThisFrame(sisIdx(1),:) - coordsThisFrame(sisIdx(2),:));
                vecN = (coordsThisFrame(sisIdx(1),1) - coordsThisFrame(sisIdx(2),1)) ./ dis;
                sisDis(idx) = dis;
                sisAngles(idx) = acos(abs(vecN(1)));
            end
            
            % sort unaligned pairs
            sisToFit = [sisToFit(~unalignedSis); sisToFit(unalignedSis)];
            centreCoords1 = [centreCoords1(~unalignedSis,:); centreCoords1(unalignedSis,:)];
            sisDis = [sisDis(~unalignedSis); sisDis(unalignedSis)];
            sisAngles = [sisAngles(~unalignedSis); sisAngles(unalignedSis)];
            nSisAligned = sum(~unalignedSis);
            
            % make cost mat with differnt cutoffs for aligned and unaligned pairs
            costMatLT = createSparseDistanceMatrix(centreCoords1(1:nSisAligned,:),centreCoords2(1:nPairsAligned,:),centreSR(1));
            costMatLB = createSparseDistanceMatrix(centreCoords1(nSisAligned+1:end,:),centreCoords2(1:nPairsAligned,:),centreSR(2));
            costMatRT = createSparseDistanceMatrix(centreCoords1(1:nSisAligned,:),centreCoords2(nPairsAligned+1:end,:),centreSR(2));
            costMatRB = createSparseDistanceMatrix(centreCoords1(nSisAligned+1:end,:),centreCoords2(nPairsAligned+1:end,:),centreSR(2));
            costMat = [costMatLT costMatRT; costMatLB costMatRB];
            costMat = full(costMat);
            costMat(costMat==0) = -1;
            
            % if no links possible continue to next iteration
            if all(costMat(:) == -1)
                progressText((nFrames-iFrame+1)/(nFrames-1),'backward processing double missing timepoints')
                continue
            end
            
            % update cost
            [idxRow,idxCol] = find(costMat ~= -1);
            if size(idxRow,1) == 1
                idxRow = idxRow';
            end
            if size(idxCol,1) == 1
                idxCol = idxCol';
            end
            extraCost = abs(1 - (sisDis(idxRow) ./ newSisDis(idxCol)));
            constraint1 = abs(sisDis(idxRow) - newSisDis(idxCol)) < 0.4;
            constraint2 = abs(sisAngles(idxRow) - newAngles(idxCol)) < 20*pi/180;
            for idx = 1:length(idxRow)
                row = idxRow(idx);
                col = idxCol(idx);
                costMat(row,col) = costMat(row,col) * extraCost(idx) * constraint1(idx) * constraint2(idx);
            end
            costMat(costMat==0) = -1;
            
            % if no links possible continue to next iteration
            if all(costMat(:) == -1)
                progressText((nFrames-iFrame+1)/(nFrames-1),'backward processing double missing timepoints')
                continue
            end
            
            % perform lap
            links = lap(costMat,[],[],1);
            links = double(links(1:size(centreCoords1,1)));
            links(links > size(centreCoords2,1)) = NaN;
            
            % get conections
            for idx = 1:size(sisToFit,1)
                iSis = sisToFit(idx);
                if ~isnan(links(idx))
                    pair = pairs(links(idx),:);
                    pair = toConsider(pair)';
                    % pair up sisters by lap
                    idxThisFrame = indexMat((iSis-1)*2+1:(iSis-1)*2+2,iFrame);
                    sisCoordsThisFrame = coordsThisFrame(idxThisFrame,:);
                    sisCoordTargetFrame = coordsTargetFrame(pair,:);
                    disMat = createDistanceMatrix(sisCoordsThisFrame,sisCoordTargetFrame);
                    pairLink = double(lap(disMat));
                    pair = pair(pairLink);
                    % fill in matrix
                    indexMat((iSis-1)*2+1:(iSis-1)*2+2,targetFrame) = pair;
                end
            end
            progressText((nFrames-iFrame+1)/(nFrames-1),'backward processing double missing timepoints')
        end
        
        %% STEP 3: LOOK FOR SISTERS THAT CAN BE MERGED
        
        % join sisters by making new tracks
        %         disp('')
        %         progressText(0)
        %         newTracks = [];
        %         for iFrame = 2:nFrames-1
        %             % get unique spot ids in this frame
        %             uniqueSpotIds = unique(indexMat(:,iFrame))';
        %             uniqueSpotIds = uniqueSpotIds(~isnan(uniqueSpotIds));
        %             for idx = uniqueSpotIds
        %                 logicIdx = indexMat(:,iFrame) == idx;
        %                 if sum(logicIdx) > 1
        %                     % make new track
        %                     newTrackIndexes = find(logicIdx);
        %                     for iLoop = 1:length(newTrackIndexes)-1
        %                         for jLoop = iLoop+1:length(newTrackIndexes)
        %                             iIdx = newTrackIndexes(iLoop);
        %                             jIdx = newTrackIndexes(jLoop);
        %                             if mod(iIdx,2) == 0
        %                                 iOtherIdx = iIdx - 1;
        %                             else
        %                                 iOtherIdx = iIdx + 1;
        %                             end
        %                             if mod(jIdx,2) == 0
        %                                 jOtherIdx = jIdx - 1;
        %                             else
        %                                 jOtherIdx = jIdx + 1;
        %                             end
        %                             newTracks = [newTracks; indexMat(iIdx,1:iFrame) indexMat(jIdx,iFrame+1:end); indexMat(iOtherIdx,1:iFrame) indexMat(jOtherIdx,iFrame+1:end); indexMat(jIdx,1:iFrame) indexMat(iIdx,iFrame+1:end); indexMat(jOtherIdx,1:iFrame) indexMat(iOtherIdx,iFrame+1:end)];%#ok<AGROW>
        %                         end
        %                     end
        %                 end
        %             end
        %             progressText((iFrame-1)/(nFrames-2),'merging tracks')
        %         end
        %         % add to current tracks
        %         indexMat = [indexMat; newTracks];%#ok<AGROW>
        %         nSis = size(indexMat,1)/2;
        
        % new methods
        % make conflict matrix
        conflicts = zeros(nSis);
        for iFrame = 2:nFrames-1
            uniqueSpotIds = unique(indexMat(:,iFrame))';
            uniqueSpotIds = uniqueSpotIds(~isnan(uniqueSpotIds));
            for idx = uniqueSpotIds
                logicIdx = indexMat(:,iFrame) == idx;
                if sum(logicIdx) > 1
                    confs = find(logicIdx);
                    confs = unique(ceil(confs/2));
                    for iLoop = 1:length(confs)-1
                        for jLoop = iLoop+1:length(confs)
                            conflicts(confs(iLoop),confs(jLoop)) = 1;
                            conflicts(confs(jLoop),confs(iLoop)) = 1;
                        end
                    end
                end
            end
        end
        
        if any(conflicts(:) == 1)
            % get disconected sub-graphs
            subGraphs = floodFillGraph(conflicts);
            idx = 0;
            while idx < length(subGraphs)
                idx = idx + 1;
                if length(subGraphs(idx).graph) == 1
                    subGraphs(idx) = [];
                    idx = idx - 1;
                end
            end
            
            for iConflict = 1:length(subGraphs)
                sistersInConflict = subGraphs(iConflict).graph;
                timesInConflict = [];
                conflictMap = {};
                for iFrame = 2:nFrames-1
                    idx = indexMat(:,iFrame);
                    keepIdx = false(length(idx),1);
                    for sis = sistersInConflict
                        keepIdx((sis-1)*2+1:(sis-1)*2+2) = true;
                    end
                    idx(~keepIdx) = NaN;
                    uniqueSpotIds = unique(idx)';
                    uniqueSpotIds = uniqueSpotIds(~isnan(uniqueSpotIds));
                    for idx = uniqueSpotIds
                        logicIdx = indexMat(:,iFrame) == idx;
                        if sum(logicIdx) > 1
                            confs = find(logicIdx);
                            toRem = [];
                            for iLoop = 1:length(confs)-1
                                for jLoop = iLoop+1:length(confs)
                                    if indexMat(confs(iLoop),iFrame+1) == indexMat(confs(jLoop),iFrame+1)
                                        toRem = [toRem; iLoop];%#ok<AGROW>
                                    end
                                end
                            end
                            toRem = unique(toRem);
                            confs(toRem) = [];
                            if length(confs) < 2
                                continue
                            end
                            confs = unique(ceil(confs/2));
                            timesInConflict = [timesInConflict iFrame];%#ok<AGROW>
                            conflictMap = [conflictMap; {confs'}];%#ok<AGROW>
                        end
                    end
                end
                if isempty(conflictMap)
                    continue
                end
                % paths
                paths = findAllPaths(conflictMap);
                timesInConflict = [timesInConflict nFrames];%#ok<AGROW>
                for i = 1:length(paths)
                    path = paths{i};
                    t = indexMat((path(1)-1)*2+1:(path(1)-1)*2+2,1:timesInConflict(1));
                    idx = 1;
                    for sis = path(2:end)
                        idx = idx + 1;
                        time = timesInConflict(idx-1:idx);
                        id = indexMat((sis-1)*2+1:(sis-1)*2+2,time(1):time(2));
                        if id(1,1) == t(1,end) || id(2,1) == t(2,end)
                            t = [t id(:,2:end)];%#ok<AGROW>
                        else
                            t = [t id([2 1],2:end)];%#ok<AGROW>
                        end
                    end
                    indexMat = [indexMat; t];%#ok<AGROW>
                end
            end
        end
        nSis = size(indexMat,1)/2;
        
        %% STEP 4: REMOVE CONFLICTS
        
        % make conflict matrix
        conflicts = zeros(nSis);
        for iFrame = 1:nFrames
            uniqueSpotIds = unique(indexMat(:,iFrame))';
            uniqueSpotIds = uniqueSpotIds(~isnan(uniqueSpotIds));
            for idx = uniqueSpotIds
                logicIdx = indexMat(:,iFrame) == idx;
                if sum(logicIdx) > 1
                    confs = find(logicIdx);
                    confs = unique(ceil(confs/2));
                    for iLoop = 1:length(confs)-1
                        for jLoop = iLoop+1:length(confs)
                            conflicts(confs(iLoop),confs(jLoop)) = 1;
                            conflicts(confs(jLoop),confs(iLoop)) = 1;
                        end
                    end
                end
            end
        end
        
        % loop until no conflicts
        while any(conflicts(:) == 1)
            % get disconected sub-graphs
            subGraphs = floodFillGraph(conflicts);
            idx = 0;
            while idx < length(subGraphs)
                idx = idx + 1;
                if length(subGraphs(idx).graph) == 1
                    subGraphs(idx) = [];
                    idx = idx - 1;
                end
            end
            
            % for each confilict island remove the sister pair with the most NaNs, if
            % there are two with hte same number then remove the sister pair with the
            % largest variance in inter-sister distance
            rowsToRem = zeros(2*length(subGraphs),1);
            for iGraph = 1:length(subGraphs)
                sis = subGraphs(iGraph).graph;
                numNaN = [];
                for iSis = sis
                    sisTracks = indexMat((iSis-1)*2+1:(iSis-1)*2+2,:);
                    sisTracks = sisTracks(:);
                    numNaN = [numNaN; sum(isnan(sisTracks))];%#ok<AGROW>
                end
                sis = sis(numNaN == max(numNaN));
                if length(sis) > 1
                    vars = [];
                    for iSis = sis
                        sisTracks = indexMat((iSis-1)*2+1:(iSis-1)*2+2,:);
                        coords = NaN(nFrames,6);
                        for t = 1:nFrames
                            if ~isempty(frameAlignment(t).alignedCoord)
                                c = frameAlignment(t).alignedCoord(:,1:3);
                                if ~isnan(sisTracks(1,t))
                                    coords(t,1:3) = c(sisTracks(1,t),:);
                                end
                                if ~isnan(sisTracks(2,t))
                                    coords(t,4:6) = c(sisTracks(2,t),:);
                                end
                            end
                        end
                        dis = normList(coords(:,1:3)-coords(:,4:6));
                        vars = [vars; nanvar(dis)];%#ok<AGROW>
                    end
                    vars(isnan(vars)) = Inf;%#ok<AGROW>
                    sis = sis(vars == max(vars));
                end
                sis = sis(1);
                rowsToRem((iGraph-1)*2+1:(iGraph-1)*2+2) = ((sis-1)*2+1:(sis-1)*2+2)';
            end
            
            indexMat(rowsToRem,:) = [];
            nSis = size(indexMat,1)/2;
            
            % make conflict matrix
            conflicts = zeros(nSis);
            for iFrame = 1:nFrames
                uniqueSpotIds = unique(indexMat(:,iFrame))';
                uniqueSpotIds = uniqueSpotIds(~isnan(uniqueSpotIds));
                for idx = uniqueSpotIds
                    logicIdx = indexMat(:,iFrame) == idx;
                    if sum(logicIdx) > 1
                        confs = find(logicIdx);
                        confs = unique(ceil(confs/2));
                        for iLoop = 1:length(confs)-1
                            for jLoop = iLoop+1:length(confs)
                                conflicts(confs(iLoop),confs(jLoop)) = 1;
                                conflicts(confs(jLoop),confs(iLoop)) = 1;
                            end
                        end
                    end
                end
            end
        end
        
        %% STEP 5: EXTRA SPOT FITTING
        
        % struct of potential new spot locations in pixels
        newSpots(1:nFrames) = struct('newCoords',[]);
        
        % forward in time
        for iFrame = 2:nFrames
            % sisters to extrapolate positions
            singles = [];
            doubles = [];
            % go over sister tracks, either single or double new spots
            for idx = 1:nSis
                sis1Id = (idx-1)*2 + 1;
                sis2Id = sis1Id + 1;
                if ~isnan(indexMat(sis1Id,iFrame-1)) && ~isnan(indexMat(sis2Id,iFrame-1)) && ~isnan(indexMat(sis1Id,iFrame)) && isnan(indexMat(sis2Id,iFrame))
                    singles = [singles; idx 2];%#ok<AGROW>
                elseif ~isnan(indexMat(sis1Id,iFrame-1)) && ~isnan(indexMat(sis2Id,iFrame-1)) && isnan(indexMat(sis1Id,iFrame)) && ~isnan(indexMat(sis2Id,iFrame))
                    singles = [singles; idx 1];%#ok<AGROW>
                elseif ~isnan(indexMat(sis1Id,iFrame-1)) && ~isnan(indexMat(sis2Id,iFrame-1)) && isnan(indexMat(sis1Id,iFrame)) && isnan(indexMat(sis2Id,iFrame))
                    doubles = [doubles; idx];%#ok<AGROW>
                end
            end
            
            % if no spots to fit continue
            if isempty(singles) && isempty(doubles)
                continue
            end
            
            % skip this frame if no spots in previous frame
            if isempty(frameAlignment(iFrame-1).alignedCoord)
                continue
            end
            
            % get current coords in this frame and previous frame, also get coord
            % system and centre of mass
            if ~isempty(frameAlignment(iFrame).alignedCoord)
                coords = frameAlignment(iFrame).alignedCoord(:,1:3);
            else
                coords = [];
            end
            previousCoords = frameAlignment(iFrame-1).alignedCoord(:,1:3);
            
            % singles
            if ~isempty(singles)
                c = zeros(size(singles,1),3);
                for id = 1:size(singles,1)
                    iSis = singles(id,1);
                    % make displacement vector
                    if singles(id,2) == 1
                        colIdx = (iSis-1)*2 + 2;
                        otherIdx = (iSis-1)*2 + 1;
                        disVec = coords(indexMat(colIdx,iFrame),:) - previousCoords(indexMat(colIdx,iFrame-1),:);
                    else
                        colIdx = (iSis-1)*2 + 1;
                        otherIdx = (iSis-1)*2 + 2;
                        disVec = coords(indexMat(colIdx,iFrame),:) - previousCoords(indexMat(colIdx,iFrame-1),:);
                    end
                    % reverse x direction if this is an anaphase frame
                    if ismember(iFrame,timesAfterAPhase)
                        disVec(1) = -disVec(1);
                    end
                    % make expected coords
                    c(id,:) = previousCoords(indexMat(otherIdx,iFrame-1),:) + disVec;
                end
                newSpots(iFrame).newCoords = [newSpots(iFrame).newCoords; c];%#ok<AGROW>
            end
            
            % doubles
            if ~isempty(doubles)
                c = zeros(2*size(doubles,1),3);
                for id = 1:size(doubles,1)
                    iSis = doubles(id);
                    % make expected coords
                    c((id-1)*2+1:(id-1)*2+2,:) = previousCoords(indexMat((iSis-1)*2+1:(iSis-1)*2+2,iFrame-1),:);
                end
                newSpots(iFrame).newCoords = [newSpots(iFrame).newCoords; c];%#ok<AGROW>
            end
        end
        
        % backward in time
        for iFrame = nFrames-1:-1:1
            % sisters to extrapolate positions
            singles = [];
            doubles = [];
            % go over sister tracks, either single or double new spots
            for idx = 1:nSis
                sis1Id = (idx-1)*2 + 1;
                sis2Id = sis1Id + 1;
                if ~isnan(indexMat(sis1Id,iFrame+1)) && ~isnan(indexMat(sis2Id,iFrame+1)) && ~isnan(indexMat(sis1Id,iFrame)) && isnan(indexMat(sis2Id,iFrame))
                    singles = [singles; idx 2];%#ok<AGROW>
                elseif ~isnan(indexMat(sis1Id,iFrame+1)) && ~isnan(indexMat(sis2Id,iFrame+1)) && isnan(indexMat(sis1Id,iFrame)) && ~isnan(indexMat(sis2Id,iFrame))
                    singles = [singles; idx 1];%#ok<AGROW>
                elseif ~isnan(indexMat(sis1Id,iFrame+1)) && ~isnan(indexMat(sis2Id,iFrame+1)) && isnan(indexMat(sis1Id,iFrame)) && isnan(indexMat(sis2Id,iFrame))
                    doubles = [doubles; idx];%#ok<AGROW>
                end
            end
            
            % if no spots to fit continue
            if isempty(singles) && isempty(doubles)
                continue
            end
            
            % skip this frame if no spots in previous frame
            if isempty(frameAlignment(iFrame+1).alignedCoord)
                continue
            end
            
            % get current coords in this frame and previous frame, also get coord
            % system and centre of mass
            if ~isempty(frameAlignment(iFrame).alignedCoord)
                coords = frameAlignment(iFrame).alignedCoord(:,1:3);
            else
                coords = [];
            end
            previousCoords = frameAlignment(iFrame+1).alignedCoord(:,1:3);
            
            % singles
            if ~isempty(singles)
                c = zeros(size(singles,1),3);
                for id = 1:size(singles,1)
                    iSis = singles(id,1);
                    % make displacement vector
                    if singles(id,2) == 1
                        colIdx = (iSis-1)*2 + 2;
                        otherIdx = (iSis-1)*2 + 1;
                        disVec = coords(indexMat(colIdx,iFrame),:) - previousCoords(indexMat(colIdx,iFrame+1),:);
                    else
                        colIdx = (iSis-1)*2 + 1;
                        otherIdx = (iSis-1)*2 + 2;
                        disVec = coords(indexMat(colIdx,iFrame),:) - previousCoords(indexMat(colIdx,iFrame+1),:);
                    end
                    % reverse x direction if this is an anaphase frame
                    if ismember(iFrame,timesAfterAPhase)
                        disVec(1) = -disVec(1);
                    end
                    % make expected coords
                    c(id,:) = previousCoords(indexMat(otherIdx,iFrame+1),:) + disVec;
                end
                newSpots(iFrame).newCoords = [newSpots(iFrame).newCoords; c];
            end
            
            % doubles
            if ~isempty(doubles)
                c = zeros(2*size(doubles,1),3);
                for id = 1:size(doubles,1)
                    iSis = doubles(id);
                    % make expected coords
                    c((id-1)*2+1:(id-1)*2+2,:) = previousCoords(indexMat((iSis-1)*2+1:(iSis-1)*2+2,iFrame+1),:);
                end
                newSpots(iFrame).newCoords = [newSpots(iFrame).newCoords; c];
            end
        end
        
        % check new coords for coords that are too close and transform to
        % pixels
        for iFrame = 1:nFrames
            c = newSpots(iFrame).newCoords;
            if isempty(c)
                continue
            end
            if ~isempty(frameAlignment(iFrame).alignedCoord)
                coords = frameAlignment(iFrame).alignedCoord(:,1:3);
                disMat = createDistanceMatrix(c,coords);
                [iRow,~] = find(disMat < 0.2);
                iRow = unique(iRow);
                c(iRow,:) = [];
            end
            if isempty(c)
                newSpots(iFrame).newCoords = c;%#ok<AGROW>
                continue
            end
            % check against previous potential new spot locations
            coords = previousNewSpots(iFrame).coords;
            if ~isempty(coords)
                disMat = createDistanceMatrix(c,coords);
                [iRow,~] = find(disMat < 0.05);
                iRow = unique(iRow);
                c(iRow,:) = [];
                if isempty(c)
                    newSpots(iFrame).newCoords = c;%#ok<AGROW>
                    continue
                end
            end
            % add spots to list
            previousNewSpots(iFrame).coords = [previousNewSpots(iFrame).coords; c];
            disMat = createSparseDistanceMatrix(c,c,0.2,0);
            [iRow,iCol] = find(disMat > 0);
            while ~isempty(iRow) && ~isempty(iCol)
                % remove the spot with the most conflicts, or if there is
                % not spot with the most conflicts remove of the two spots
                % with the smallest distance the spot with the smallest
                % average distance to other spots
                tmp = [iRow; iCol];
                uniqueTmp = unique(tmp);
                countTmp = zeros(1,length(uniqueTmp));
                for iTmp = 1:length(uniqueTmp)
                    countTmp(iTmp) = sum(uniqueTmp(iTmp)==tmp);
                end
                maxTmp = max(countTmp);
                numMaxTmp = sum(maxTmp==countTmp);
                if numMaxTmp == 1
                    id = maxTmp==countTmp;
                    id = uniqueTmp(id);
                    c(id,:) = [];
                else
                    disMat(disMat==0) = Inf;
                    [iRow,iCol] = find(disMat == min(disMat(:)),1);
                    disMat = createDistanceMatrix(c,c);
                    meanDisI = mean(disMat(iRow,:));
                    meanDisJ = mean(disMat(iCol,:));
                    if meanDisI < meanDisJ
                        c(iRow,:) = [];
                    else
                        c(iCol,:) = [];
                    end
                end
                disMat = createSparseDistanceMatrix(c,c,0.2,0);
                [iRow,iCol] = find(disMat > 0);
            end
            coordSystem = frameAlignment(iFrame).coordSystem;
            com = frameAlignment(iFrame).centerOfMass;
            t = 0;
            loop = 1;
            while isnan(com(1)) && loop
                t = t + 1;
                if iFrame-t < 1 && iFrame+t > nFrames
                    loop = 0;
                end
                if iFrame-t > 0
                    com = frameAlignment(iFrame-t).centerOfMass;
                    coordSystem = frameAlignment(iFrame-t).coordSystem;
                end
                if isnan(com(1))
                    if iFrame+t <= nFrames
                        com = frameAlignment(iFrame+t).centerOfMass;
                        coordSystem = frameAlignment(iFrame+t).coordSystem;
                    end
                end
            end
            frameAlignment(iFrame).centerOfMass = com;
            frameAlignment(iFrame).coordSystem = coordSystem;
            % transform to raw coords and transform into pixles
            c = (coordSystem * c')';
            c = c + repmat(com,size(c,1),1);
            c = c ./ repmat(pixelSize,size(c,1),1);
            % make sure coords aren't outsize the frame
            c(c(:,1)<=0,1) = 0.0001;
            c(c(:,2)<=0,2) = 0.0001;
            c(c(:,3)<=0,3) = 0.0001;
            c(c(:,1)>xSize,1) = xSize;
            c(c(:,2)>ySize,2) = ySize;
            c(c(:,3)>zSize,3) = zSize;
            newSpots(iFrame).newCoords = c;%#ok<AGROW>
        end
        
        
        % perform fitting, change the order of spots afterwards to match current
        % indexes, flag whether or not any new spot have been ades at all
        disp('')
        progressText(0,'re-fitting for new spots: xxx new spots found in frame xxx')
        newSpotAdded = false;
        for iFrame = 1:nFrames
            newSpotsThisFrame = false;
            
            % if no new coords continue
            if isempty(newSpots(iFrame).newCoords)
                progressText(iFrame/nFrames,['re-fitting for new spots: no new spots for frame ' int2str(iFrame)])
                continue
            end
            if isempty(initCoord(iFrame).allCoordPix)
                coords = newSpots(iFrame).newCoords;
            else
                coords = [initCoord(iFrame).allCoordPix(:,1:3); newSpots(iFrame).newCoords];
            end
            % fit
            [coordList,ampList,bgList] = spotMMFit(movie(:,:,:,iFrame),coords,psfSigma,'influenceRadius',3*psfSigma,'overlapFactor',1000*psfSigma);
            
            if isempty(coordList)
                progressText(iFrame/nFrames,['re-fitting for new spots: no new spots for frame ' int2str(iFrame)])
                continue
            end
            
            if size(coordList,1) ~= initCoord(iFrame).nSpots
                newSpotsThisFrame = true;
            end
            
            if ~isempty(initCoord(iFrame).allCoordPix)
                disMat = createSparseDistanceMatrix(initCoord(iFrame).allCoordPix(:,1:3),coordList(:,1:3),1);
                links = lap(disMat,[],[],1);
                links = double(links(1:initCoord(iFrame).nSpots));
                % remove any previously removed spots
                while any(links > size(coordList,1))
                    newSpotsThisFrame = true;
                    idx = find(links > size(coordList,1),1);
                    % remove this idx from indexMat
                    indexMat(indexMat(:,iFrame) == idx,iFrame) = NaN;
                    % decrese the index of other spots with higher indexes
                    indexMat(indexMat(:,iFrame) > idx,iFrame) = indexMat(indexMat(:,iFrame) > idx,iFrame) - 1;
                    % remove from links
                    links(idx) = [];
                end
                otherSpots = setxor(1:size(coordList,1),links);
                coordList = [coordList(links,:); coordList(otherSpots,:)];
                ampList = [ampList(links,:); ampList(otherSpots,:)];
                bgList = [bgList(links,:); bgList(otherSpots,:)];
            else
                newSpotsThisFrame = true;
            end
            
            if ~newSpotsThisFrame
                progressText(iFrame/nFrames,['re-fitting for new spots: no new spots for frame ' int2str(iFrame)])
                continue
            end
            
            newSpotAdded = true;
            
            progressText(iFrame/nFrames,['re-fitting for new spots: ' int2str(size(coordList,1)-initCoord(iFrame).nSpots) ' new spots for frame ' int2str(iFrame)])
            
            % add back into initCoord
            initCoord(iFrame).nSpots = size(coordList,1);
            initCoord(iFrame).allCoordPix = coordList;
            initCoord(iFrame).amp = ampList;
            initCoord(iFrame).bg = bgList;
            initCoord(iFrame).allCoord = initCoord(iFrame).allCoordPix .* repmat(pixelSize,initCoord(iFrame).nSpots,2);
        end
        
        % if no spots added end loop
        if ~newSpotAdded
            fittingLoop1 = false;
        else
            
            % do a plane fit to detect outliers
            dataStruct_tmp = dataStruct;
            dataStruct_tmp.initCoord = initCoord;
            if ~isempty(planeFit)
                dataStruct_tmp = makiFitPlane(dataStruct_tmp,0);
            end
            
            % update outlier lists and update coords
            for iFrame = 1:nFrames
                if initCoord(iFrame).nSpots > 0
                    if ~isempty(planeFit)
                        updatedClass(iFrame).unalignedIdx = dataStruct_tmp.planeFit(iFrame).unalignedIdx;
                    end
                    coordSystem = frameAlignment(iFrame).coordSystem;
                    com = frameAlignment(iFrame).centerOfMass;
                    coords = initCoord(iFrame).allCoord;
                    coords(:,1:3) = coords(:,1:3) - repmat(com,size(coords,1),1);
                    coords(:,1:3) = (coordSystem \ coords(:,1:3)')';
                    frameAlignment(iFrame).alignedCoord = coords;
                end
            end
            
            
        end
    end
    
    %% STEP 6: track remaining spots
    
    % remove spots that have been assigned to sisters, save to a new initCoord,
    % fameAlignment and updatedClass
    initCoord_tmp = initCoord;
    frameAlignment_tmp = frameAlignment;
    updatedClass_tmp = updatedClass;
    % store the indexes of spots kept
    keptSpots(1:nFrames) = struct('idx',[]);
    for iFrame = 1:nFrames
        if initCoord_tmp(iFrame).nSpots > 0
            spotIdx = indexMat(:,iFrame);
            spotIdx = spotIdx(~isnan(spotIdx));
            coords = frameAlignment(iFrame).alignedCoord;
            keepIdx = setxor(1:size(coords,1),spotIdx);
            keptSpots(iFrame).idx = keepIdx;%#ok<AGROW>
            unalignedIdx = updatedClass(iFrame).unalignedIdx;
            laggingIdx = updatedClass(iFrame).laggingIdx;
            unalignedIdx = find(ismember(keepIdx,unalignedIdx));
            laggingIdx = find(ismember(keepIdx,laggingIdx));
            alignedIdx = setxor(1:length(keepIdx),[unalignedIdx laggingIdx]);
            coords = coords(keepIdx,:);
            frameAlignment_tmp(iFrame).alignedCoord = coords;
            updatedClass_tmp(iFrame).unalignedIdx = unalignedIdx;
            updatedClass_tmp(iFrame).laggingIdx = laggingIdx;
            updatedClass_tmp(iFrame).inlierIdx = alignedIdx;
            initCoord_tmp(iFrame).allCoord = initCoord_tmp(iFrame).allCoord(keepIdx,:);
            if ~isempty(planeFit)
                planeFit(iFrame).rotatedCoord = frameAlignment_tmp(iFrame).alignedCoord;
            end
        end
    end
    
    % gereate tracks from remaining spots
    dataStruct_tmp = dataStruct;
    dataStruct_tmp.frameAlignment = frameAlignment_tmp;
    dataStruct_tmp.updatedClass = updatedClass_tmp;
    dataStruct_tmp.initCoord = initCoord_tmp;
    if ~isempty(planeFit)
        dataStruct_tmp.planeFit = planeFit;
    end
    disp('re-tracking remaining spots')
    dataStruct_tmp = makiGenerateTracks_frameAlignment(dataStruct_tmp);
    
    % update phase
    for iFrame = timesAfterAPhase
        dataStruct_tmp.planeFit(iFrame).phase = 'a';
    end
    
    % group sisters
    disp('checking for any more sisters')
    dataStruct_tmp = makiGroupSisters(dataStruct_tmp);
    
    % if there are any new sisters add the right index to the matrix and begin
    % all over again...
    if length(dataStruct_tmp.sisterList) == 1 && isempty(dataStruct_tmp.sisterList.distances)
        fittingLoop2 = false;
        disp('no new sisters found, stopping...')
    else
        % else add to matrix
        sList = dataStruct_tmp.sisterList(1).trackPairs(:,1:2);
        
        extraSis = [];
        for iSis = 1:size(sList,1)
            
            % get coords and amps
            coords1 = getCoords(dataStruct_tmp.tracks(sList(iSis,1)),nFrames);
            coords1(:,4) = getAmp(dataStruct_tmp.tracks(sList(iSis,1)),nFrames);
            coords2 = getCoords(dataStruct_tmp.tracks(sList(iSis,2)),nFrames);
            coords2(:,4) = getAmp(dataStruct_tmp.tracks(sList(iSis,2)),nFrames);
            
            % detect outliers
            outliers1 = 0;
            outliers2 = 0;
            while ~isempty(outliers1) || ~isempty(outliers2)
                [outliers1,outliers2] = detectOuliers(coords1,coords2,timesBeforeAphase,maxSisDis);
                
                % remove outliers (make NaN)
                coords1(outliers1,:) = NaN;
                coords2(outliers2,:) = NaN;
            end
            % find not NaNs on either sister
            coordIdx = ~isnan(coords1(:,1)) | ~isnan(coords2(:,1));
            
            % cluster coords
            clusters = bwconncomp([coordIdx coordIdx]);
            
            % get original indexes
            featIdx1 = getFeatIdx(dataStruct_tmp.tracks(sList(iSis,1)),nFrames)';
            featIdx2 = getFeatIdx(dataStruct_tmp.tracks(sList(iSis,2)),nFrames)';
            
            % convert to proper spot idxs
            for iFrame = 1:nFrames
                if ~isnan(featIdx1(iFrame))
                    idx = keptSpots(iFrame).idx;
                    featIdx1(iFrame) = idx(featIdx1(iFrame));
                end
                if ~isnan(featIdx2(iFrame))
                    idx = keptSpots(iFrame).idx;
                    featIdx2(iFrame) = idx(featIdx2(iFrame));
                end
            end
            
            % cluster indicies and add to matrix
            for iCluster = 1:clusters.NumObjects
                idx = clusters.PixelIdxList{iCluster}(1:end/2);
                idx = setxor(1:nFrames,idx);
                featIdx1_tmp = featIdx1;
                featIdx1_tmp(idx) = NaN;
                featIdx2_tmp = featIdx2;
                featIdx2_tmp(idx) = NaN;
                extraSis = [extraSis; featIdx1_tmp; featIdx2_tmp];%#ok<AGROW>
            end
        end
        
        % check againt previous extra sisters
        for iExSis = 1:length(previousExtraSis)
            pExtraSis = previousExtraSis(iExSis).pExtraSis;
            if all(size(pExtraSis) == size(extraSis))
                previousExtraSis_tmp = pExtraSis;
                extraSis_tmp = extraSis;
                previousExtraSis_tmp(isnan(previousExtraSis_tmp)) = 0;
                extraSis_tmp(isnan(extraSis_tmp)) = 0;
                if all(previousExtraSis_tmp(:) == extraSis_tmp(:))
                    disp('same sisters found, stopping...')
                    fittingLoop2 = false;
                    break
                end
            end
        end
        clear previousExtraSis_tmp
        previousExtraSis_tmp.pExtraSis = extraSis;
        previousExtraSis = [previousExtraSis previousExtraSis_tmp];%#ok<AGROW>
        if fittingLoop2
            disp('more sisters found, continuing...')
            indexMat = [indexMat; extraSis];%#ok<AGROW>
        end
    end
    nSis = size(indexMat,1)/2;
end

%% STEP 7: ASSEMBLE NEW TRACKS AND NEW SISTERLIST

% check and remove any all NaN sisters
toRem = [];
for idx = 1:size(indexMat,1)
    if all(isnan(indexMat(idx,:)))
        if mod(idx,2) == 0
            otherIdx = idx - 1;
        else
            otherIdx = idx + 1;
        end
        toRem = [toRem; idx; otherIdx];%#ok<AGROW>
    end
end
toRem = unique(toRem);
indexMat(toRem,:) = [];

nSis = size(indexMat,1)/2;

% add extra tracks to the end of the matrix
for iTrack = 1:length(dataStruct_tmp.tracks)
    featIdx = getFeatIdx(dataStruct_tmp.tracks(iTrack),nFrames);
    for iFrame = 1:length(featIdx)
        if ~isnan(featIdx(iFrame))
            idx = keptSpots(iFrame).idx;
            featIdx(iFrame) = idx(featIdx(iFrame));
        end
    end
    indexMat = [indexMat; featIdx'];%#ok<AGROW>
end

% make tracks
clear tracks
tracks(1:size(indexMat,1)) = struct('tracksFeatIndxCG',[],'tracksCoordAmpCG',[],'seqOfEvents',[1 1 1 NaN; 1 2 1 NaN],'coordAmp4Tracking',[],'dependencies',[]);
for iTrack = 1:size(indexMat,1)
    index = indexMat(iTrack,:);
    startTime = find(~isnan(index),1);
    endTime = find(~isnan(index),1,'last');
    times = startTime:endTime;
    tracksCoordAmpCG = NaN(1,8*length(times));
    coordAmp4Tracking = NaN(1,8*length(times));
    index = index(startTime:endTime);
    index(isnan(index)) = 0;
    tracks(iTrack).tracksFeatIndxCG = index;
    tracks(iTrack).seqOfEvents(1,1) = startTime;
    tracks(iTrack).seqOfEvents(2,1) = endTime;
    count = 0;
    for iFrame = times
        count = count + 1;
        if index(count) ~= 0
            coords = initCoord(iFrame).allCoord;
            coordsA = frameAlignment(iFrame).alignedCoord;
            amp = initCoord(iFrame).amp;
            tracksCoordAmpCG((count-1)*8+1) = coords(index(count),1);
            tracksCoordAmpCG((count-1)*8+2) = coords(index(count),2);
            tracksCoordAmpCG((count-1)*8+3) = coords(index(count),3);
            tracksCoordAmpCG((count-1)*8+4) = amp(index(count),1);
            tracksCoordAmpCG((count-1)*8+5) = coords(index(count),4);
            tracksCoordAmpCG((count-1)*8+6) = coords(index(count),5);
            tracksCoordAmpCG((count-1)*8+7) = coords(index(count),6);
            tracksCoordAmpCG((count-1)*8+8) = amp(index(count),2);
            coordAmp4Tracking((count-1)*8+1) = coordsA(index(count),1);
            coordAmp4Tracking((count-1)*8+2) = coordsA(index(count),2);
            coordAmp4Tracking((count-1)*8+3) = coordsA(index(count),3);
            coordAmp4Tracking((count-1)*8+4) = amp(index(count),1);
            coordAmp4Tracking((count-1)*8+5) = coordsA(index(count),4);
            coordAmp4Tracking((count-1)*8+6) = coordsA(index(count),5);
            coordAmp4Tracking((count-1)*8+7) = coordsA(index(count),6);
            coordAmp4Tracking((count-1)*8+8) = amp(index(count),2);
        end
    end
    tracks(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
    tracks(iTrack).coordAmp4Tracking = coordAmp4Tracking;
end

% make sisterList
clear sisterList
sisterList(1:nSis) = struct('trackPairs',[],'coords1',[],'coords2',[],'sisterVectors',[],'distances',[],'dependencies',[]);
sisV = (1:nSis)';
sisterList(1).trackPairs = [(sisV-1)*2+1 (sisV-1)*2+2];
for iSis = 1:nSis
    coords1 = NaN(nFrames,6);
    coords2 = NaN(nFrames,6);
    times = tracks(sisterList(1).trackPairs(iSis,1)).seqOfEvents(1,1):tracks(sisterList(1).trackPairs(iSis,1)).seqOfEvents(2,1);
    count = 0;
    for iFrame = times
        count = count + 1;
        coords1(iFrame,1:3) = tracks(sisterList(1).trackPairs(iSis,1)).coordAmp4Tracking((count-1)*8+1:(count-1)*8+3);
        coords1(iFrame,4:6) = tracks(sisterList(1).trackPairs(iSis,1)).coordAmp4Tracking((count-1)*8+5:(count-1)*8+7);
    end
    times = tracks(sisterList(1).trackPairs(iSis,2)).seqOfEvents(1,1):tracks(sisterList(1).trackPairs(iSis,2)).seqOfEvents(2,1);
    count = 0;
    for iFrame = times
        count = count + 1;
        coords2(iFrame,1:3) = tracks(sisterList(1).trackPairs(iSis,2)).coordAmp4Tracking((count-1)*8+1:(count-1)*8+3);
        coords2(iFrame,4:6) = tracks(sisterList(1).trackPairs(iSis,2)).coordAmp4Tracking((count-1)*8+5:(count-1)*8+7);
    end
    sisterVectors = coords1(:,1:3) - coords2(:,1:3);
    sisterVectors(:,4:6) = sqrt(coords1(:,4:6).^2 + coords2(:,4:6).^2);
    distances = normList(sisterVectors(:,1:3));
    distances(:,2) = sqrt(sum((sisterVectors(:,1:3)./repmat(distances,1,3)).^2 .* sisterVectors(:,4:6).^2,2));
    sisterList(iSis).coords1 = coords1;
    sisterList(iSis).coords2 = coords2;
    sisterList(iSis).sisterVectors = sisterVectors;
    sisterList(iSis).distances = distances;
end

% add to dataStruct
dataStruct.initCoord = initCoord;
dataStruct.tracks = tracks;
dataStruct.sisterList = sisterList;

% re-fitting plane
if ~isempty(planeFit)
    disp('re-fiting plane')
    dataStruct = makiFitPlane(dataStruct,0);
end
% re-updating classes
disp('re-updating classes')
dataStruct = makiUpdateClass(dataStruct);
% re-align frames
disp('re-aligning frames')
dataStruct = makiAlignAllFrames(dataStruct);


%% improve anaphase detection (again)
% improve anaphase detection by looking for deviations in average sister
% separation



% get aligned sister coords
sisterList = makiConstructAlignedSisters_allCoords(dataStruct);

% get cell phase
phase = catStruct(1,'dataStruct.updatedClass.phase',[],1); % phases of the cell cycle
firstFrameAPhase = find(phase=='a',1); % find the first frame theat the cell goes into aphase, if it does at all
if isempty(firstFrameAPhase) % if the cell did not go into aphase
    endMetaphase = nFrames; % test tracks until the end of the movie
else
    endMetaphase = firstFrameAPhase - 1; % else test upto the time just before aphase
end

% only if there was a planeFit
if ~isempty(planeFit)
    
    % conpair to previous gess of the start of metaphase
    endMetaphase = min(endMetaphase,max(timesBeforeAphase));
    
    avSep = zeros(endMetaphase,1);
    for t = 1:endMetaphase
        tmp = NaN(nSis,1);
        for iSis = 1:nSis
            coords1 = sisterList(iSis).coords1(t,1:3);
            coords2 = sisterList(iSis).coords2(t,1:3);
            tmp(iSis) = normList(coords1-coords2);
        end
        avSep(t) = robustMean(tmp);
    end
    
    % take mean
    avSepM = robustMean(avSep);
    
    % do robust exp fitting
    [~,~,idx] = robustExponentialFit2(avSep);
    idx = setxor(1:endMetaphase,idx);
    
    % go though outliers, if they are at the end of the movie and their sep is
    % above average then make them anaphase frames
    newA = [];
    if ~isempty(idx)
        loop = 1;
        count = -1;
        while loop && count < endMetaphase
            count = count + 1;
            if idx(end) == endMetaphase - count && avSep(idx(end)) > avSepM
                newA = [idx(end); newA];%#ok<AGROW>
                idx(end) = [];
            else
                loop = 0;
            end
        end
    end
    
else
    newA = [];
end

% get times before and after A phase
timesBeforeAphase = setxor(1:endMetaphase,newA)';
timesAfterAPhase = setxor(1:nFrames,timesBeforeAphase)';
if ~isempty(timesAfterAPhase)
    timesAfterAPhase = returnRightVector(timesAfterAPhase,[],'r');
else
    timesAfterAPhase = [];
end

% update phase
for iFrame = timesAfterAPhase
    dataStruct.planeFit(iFrame).phase = 'a';
    dataStruct.updatedClass(iFrame).phase = 'a';
end


%% SUBFUNCTIONS

    function [outliers1,outliers2] = detectOuliers(coords1, coords2, timesToTest, maxSisDis)
        
        % crop down coords
        coords1 = coords1(timesToTest,:);
        coords2 = coords2(timesToTest,:);
        
        % separation
        sep = normList(coords1 - coords2);
        
        % outlier detection using least-median squares with a cuttoff of 5
        % sigma for x and y and z, strategy: outlier detection in change in
        % coords, absolute cutoff of 1.8 microns for sep
        outliers1 = [];
        outliers2 = [];
        [~,~,~,tmp] = robustMean(diff(coords1(:,1)),[],5);
        outliers1 = [outliers1; tmp+1];
        [~,~,~,tmp] = robustMean(diff(coords1(:,2)),[],5);
        outliers1 = [outliers1; tmp+1];
        [~,~,~,tmp] = robustMean(diff(coords1(:,3)),[],5);
        outliers1 = [outliers1; tmp+1];
        % robustExp for amp
        idx = find(~isnan(coords1(:,1)));
        [~,~,tmp] = robustExponentialFit2(coords1(:,4));
        tmp = idx(tmp);
        tmp = setxor(1:size(coords1,1),tmp);
        tmp = setdiff(tmp,find(isnan(coords1(:,1))));
        if ~isempty(tmp)
            tmp = returnRightVector(tmp,[],'c');
        else
            tmp = [];
        end
        outliers1 = [outliers1; tmp];
        [~,~,~,tmp] = robustMean(diff(coords2(:,1)),[],5);
        outliers2 = [outliers2; tmp+1];
        [~,~,~,tmp] = robustMean(diff(coords2(:,2)),[],5);
        outliers2 = [outliers2; tmp+1];
        [~,~,~,tmp] = robustMean(diff(coords2(:,3)),[],5);
        outliers2 = [outliers2; tmp+1];
        % robustExp for amp
        idx = find(~isnan(coords2(:,1)));
        [~,~,tmp] = robustExponentialFit2(coords2(:,4));
        tmp = idx(tmp);
        tmp = setxor(1:size(coords2,1),tmp);
        tmp = setdiff(tmp,find(isnan(coords2(:,1))));
        if ~isempty(tmp)
            tmp = returnRightVector(tmp,[],'c');
        else
            tmp = [];
        end
        outliers2 = [outliers2; tmp];
        tmp = find(sep > maxSisDis);
        outliers1 = [outliers1; tmp];
        outliers2 = [outliers2; tmp];
        
        % keep unique results
        outliers1 = unique(outliers1);
        outliers2 = unique(outliers2);
        
        % update outliers according to testing frames
        outliers1 = timesToTest(outliers1);
        outliers2 = timesToTest(outliers2);
    end

    function featIdx = getFeatIdx(track,noTimePoints)
        % gets the featIdx of a track for all times points, leaving a NaN
        % for no feat
        
        featIdx = NaN(noTimePoints,1);
        
        startTime = track.seqOfEvents(1,1); % get track start and end times
        endTime = track.seqOfEvents(2,1);
        
        featIdx(startTime:endTime) = track.tracksFeatIndxCG';
        
        featIdx(featIdx==0) = NaN;
        
    end

    function coords = getCoords(track,noTimePoints)
        % gets the featIdx of a track for all times points, leaving a NaN
        % for no feat
        
        coords = NaN(noTimePoints,3);
        
        startTime = track.seqOfEvents(1,1); % get track start and end times
        endTime = track.seqOfEvents(2,1);
        
        coords(startTime:endTime,1) = track.coordAmp4Tracking(1:8:end)';
        coords(startTime:endTime,2) = track.coordAmp4Tracking(2:8:end)';
        coords(startTime:endTime,3) = track.coordAmp4Tracking(3:8:end)';
        
    end

    function amp = getAmp(track,noTimePoints)
        amp = NaN(noTimePoints,1);
        startTime = track.seqOfEvents(1,1); % get track start and end times
        endTime = track.seqOfEvents(2,1);
        amp(startTime:endTime) = track.tracksCoordAmpCG(4:8:end)';
    end


end

