function dataStruct = makiFitPlaneCENPQ(dataStruct,verbose)
% EHarry Sept 2011

%% ORIGINAL HEADER
% % %MAKIFITPLANE attempts to fit planes into kinetochore clusters
% % %
% % % SYNOPSIS: dataStruct = makiFitPlane(dataStruct,verbose)
% % %
% % % INPUT dataStruct (opt): maki data structure. If empty, it will be loaded
% % %                         via GUI
% % %       verbose (opt)   : 0: No plotting
% % %                         1: plot results in the end (default)
% % %                         2: also plot plane fit for every frame
% % %                         The plots of verbose==1 can also be generated by
% % %                         calling makiFitPlanePlot(dataStruct)
% % %
% % % OUTPUT dataStruct.planeFit structure of length nTimepoints. 
% % %
% % %               for a timepoint where no plane could either be derived from the
% % %               eigenvectors or interpolated only the fields .planeVectorClassifier,
% % %               .eigenVectors, .eigenValues and .planeOrigin (which is really
% % %               the center of mass) are populated. 
% % %
% % %               .plane [a b c d] for the plane equation ax+by+cz-d=0; 
% % %               .planeCoord      spot coordinates relative to the plane.
% % %                   The first dimension is perpendicular to the plane, the
% % %                   second direction is the intersect between the plane and
% % %                   the xy-plane
% % %               .planeVectorClassifier 1 | 0 dependent on whether the
% % %                   planeVectors orginate from eigenvalues, or are
% % %                   interpolated (onset of anaphase). In prophase and early
% % %                   prometaphase, where planeVectors can neither be derived
% % %                   from eigenvalues nor interpolated, the classifier is
% % %                   set to 0 and the .planeVectors field is empty
% % %               .planeVectors    unit vectors of the plane coordinate
% % %                   system
% % %               .planeOrigin     Origin of the plane (center of mass of the
% % %                   inlier spots)
% % %               .eigenVectors    eigenVectors of the spot covariance matrix
% % %               .eigenValues     corresponding eigenValues
% % %               .inlierIdx       index into initCoord of all spots that are
% % %                   considered to be inliers
% % %               .unalignedIdx    index into initCoord of all spots that are
% % %                   considered to belong to lagging chromosomes 
% % %                   (occurs late prometaphase through anaphase)
% % %               .laggingIdx      index into initCoord of all spots that are
% % %                   considered to belong to lagging chromosomes
% % %                   (occurs only in anaphase frames)
% % %               .phase 
% % %                   'e': prophase/early prometaphase -> no plane fit
% % %                        possible;
% % %                   'p': late prometaphase -> plane fit possible but 
% % %                        unaligned kinetochores found; 
% % %                   'm': metaphase -> plane fit possible and no unaligned 
% % %                        kinetochores found;
% % %                   'a': anaphase -> either planefit possible with
% % %                        eigenvalue normal > mean eigenvalue, or no
% % %                        plane fit possible but frame comes after metaphase
% % %                        or anaphase frames (happens when confusing
% % %                        anaphase frames come at the end of a movie).
% % %               .distParms       [variance,skewness,kurtosis,pNormal]' for
% % %                   the planeCoord. pNormal is the p-value that the data
% % %                   comes from a normal distribution. The tabulated values
% % %                   for the Lilliefors test only go from 0.01 to 0.5!
% % %               .deltaP          p-value of the ks-test comparing the
% % %                   distribution of the distances of inlier spots from the
% % %                   metaphase plate of t and t+1
% % %               .deltaAngle      difference in angle between planes in t
% % %                   and t+1
% % %
% % % created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
% % %
% % % created by: jdorn, kjaqaman, gdanuser
% % % DATE: 03-Jul-2007
% % %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warningState = warning;
warning off stats:lillietest:OutOfRangeP

% threshold for an acceptable eigenvalue ratio
minEigenValRatio1 = 3;
minEigenValRatio2 = 3;

% rank of n of neighbors used for the initial detection of outlier
% kinetochores; n/2 reflects the number of expected unaligned sister pairs
rankNNearestNeighbors = 40; 

% minimal number of consecutive frames in a movie that have stable enough
% eigenvectors for plane rotation estimation
minConsecFrames = 5; 

%TEST input
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end
if nargin < 2 || isempty(verbose)
    verbose = 1;
end

%determine whether to use 2D projection instead of 3D
%the 2D approximation only works in late prometaphase and metaphase
dataProperties = dataStruct.dataProperties;
if ~isfield(dataProperties,'planeFitParam')
    use2D = 0;
else
    use2D = dataProperties.planeFitParam.use2D;
end

% check whether analysis has been done
if isempty(dataStruct.initCoord)
    dataStruct = makiInitCoord(dataStruct,0);
end

% get coordinates, dataProperties, etc
initCoord = dataStruct.initCoord;
nTimePoints = length(initCoord);
nSpots = cat(1,initCoord.nSpots);

%assign dimensionality for fit
probDim = 3 - use2D;

planeFit(1:nTimePoints) = struct('plane',[],'planeCoord',[],'planeVectorClassifier', 0, ...
    'planeVectors',[],'planeOrigin',[],'eigenVectors',[],'eigenValues',[],...
    'inlierIdx',[],'unalignedIdx',[],'laggingIdx',[],'phase','e',...
    'distParms',[],'deltaP',[],'deltaAngle',[]);

% initialize lists of frames with and without plane
framesNoPlane = 1:nTimePoints;
framesWiPlane = [];

% loop through timepoints. Get covariance of point cloud, and the
% corresponding eigenvalues. Label frames that have sufficient anisotropy.

eigenValues = zeros(nTimePoints,probDim);
eigenVectors = zeros(probDim,probDim,nTimePoints);  %vectors in cols
meanCoord = zeros(nTimePoints,probDim);
meanCoordFull = zeros(nTimePoints,3);

goodFrames1 = [];
potFrames = [];
for t=1:nTimePoints

    % make an initial outlier detection
    dist = pdist(initCoord(t).allCoord(:,1:probDim));
    dist = squareform(dist); 
    sortedDist = sort(dist);
    meanNearestDist = mean(sortedDist(2:rankNNearestNeighbors+1,:));
    
    % get kinetochores whose mean distance to the rankNNearestNeighbors
    % falls within the scatter of most of the kinetochores
    [dummy, dummy, inlierIdx, outlierIdx] = robustMean(meanNearestDist,2);
    
    % store an initial inlier index (this may be modified later)
    planeFit(t).inlierIdx = inlierIdx;
    planeFit(t).unalignedIdx = outlierIdx;
    
    % tree = linkage(dist);
    % dendrogram(tree);
    % treeInc = inconsistent(tree,size(initCoord(t).allCoord,1)-1); 
    
    % get data for eigenRatio, subsequent plane fitting
    [eigenVectors(:,:,t), eigenValues(t,:), meanCoord(t,:)] = ...
        eigenCalc(initCoord(t).allCoord(inlierIdx,1:probDim));
    [dummy,dummy,meanCoordFull(t,:)] = ...
        eigenCalc(initCoord(t).allCoord(inlierIdx,1:3)); % needed only for 3D center of mass
    
    % fill in the center of mass into the planeOrigin no matter whether
    % there will ever be a plane
    planeFit(t).planeOrigin = meanCoordFull(t,:);
    planeFit(t).eigenVectors = eigenVectors(:,:,t);
    planeFit(t).eigenValues = eigenValues(t,:);
    
    % classify the anisotropy of the point cloud
    [maxEigenVal, maxIndx] = max(eigenValues(t,:));
    [minEigenVal, minIndx] = min(eigenValues(t,:));
    
    % if there is sufficient anisotropy, classify this frame as a good frame 
    if( (eigenValues(t,maxIndx)/mean(eigenValues(t,setdiff(1:probDim,maxIndx))) ...
            > minEigenValRatio1) ||...
            (eigenValues(t,minIndx)/mean(eigenValues(t,setdiff(1:probDim,minIndx))) ...
            < 1/minEigenValRatio1))
        planeFit(t).planeVectorClassifier = 1;
        goodFrames1 = [goodFrames1 t];
    end
    
    %indicate if this frame has borderline anisotropy
    if( (eigenValues(t,maxIndx)/mean(eigenValues(t,setdiff(1:probDim,maxIndx))) ...
            > minEigenValRatio2) ||...
            (eigenValues(t,minIndx)/mean(eigenValues(t,setdiff(1:probDim,minIndx))) ...
            < 1/minEigenValRatio2))
        potFrames = [potFrames t];
    end
    
end %(for t=1:nTimePoints)

%if there are enough good frames, go over potentially good frames and
%convert them into good frames if they come after good frames
goodFrames = goodFrames1;
if length(goodFrames) >= minConsecFrames
    potFrames = setdiff(potFrames,goodFrames1);
    for t = potFrames
        if any(goodFrames1<t)
            goodFrames = [goodFrames t];
            planeFit(t).planeVectorClassifier = 1;
        end
    end
    goodFrames = sort(goodFrames);
end

% we set the condition that at least minConsecFrames pairs of consecutive frames must
% exist, where the eignevectors are sufficiently stable. Otherwise, the
% movie is considered as a prophase/early prometaphase movie with no plane
nConsecFrames = length(find(diff(goodFrames) == 1));
if nConsecFrames >= minConsecFrames

    % assign the eigenvectors so that they generate minimal global rotation
    % between good frames
    [eigenVecAssign,eigenVectors(:,:,goodFrames),eigenVectorRotCos] = assignEigenVecs(eigenVectors(:,:,goodFrames));

    % assume no rotation between the first frame of the movie and the virtual
    % time point before
    eigenVectorCos = [ones(probDim,1) eigenVectorRotCos];
    eigenVectorRotation = acos(eigenVectorCos);

    if ~use2D

        evecScore= [0;0;0];

        % define the normal vector of the rotating plane as the one whose eigenvalue is
        % overall the most distant from the two other eigenvalues
        for t = 1:length(goodFrames)

            % copy updated eigenvectors and eigenvalues into the data
            % structure
            planeFit(goodFrames(t)).eigenVectors = eigenVectors(:,:,goodFrames(t));
            planeFit(goodFrames(t)).eigenValues = eigenValues(goodFrames(t),:);

            % calculate geometric means of the pairwise distances in every time
            % point
            diffs = pdist(eigenValues(goodFrames(t),:)');
            geomDist(1,t) = sqrt(diffs(1)*diffs(2));
            geomDist(2,t) = sqrt(diffs(1)*diffs(3));
            geomDist(3,t) = sqrt(diffs(2)*diffs(3));

            % score : minimize the rotation of the normal and maximize the geometric
            % mean distance of the eigenvalue associated with the normal to the two
            % in plane eignevalues.
            % on average: eigenvectors associated with the plane normal have a
            % large geometric mean difference to inplane eigenvectors, and inplane
            % eigenvectors tend to rotate more
            evecScore(:) = evecScore(:) + eigenVectorRotation(:,t)./geomDist(eigenVecAssign(:,t),t);

        end

        % the normal is the one with the largest distance
        [dummy, normalIndx] = min(evecScore);

    else

        % the above does not make sense in 2D. I will assume that the
        % cases that need 2D are always in late prometaphase or metaphase, in
        % which case the smaller eigenvalue is the normal - KJ
        evecScore = zeros(1,2);
        for t = 1 : length(goodFrames)
            evecScore = evecScore + eigenValues(goodFrames(t),eigenVecAssign(:,goodFrames(t)));
        end
        [dummy,normalIndx] = min(evecScore);

    end

    %before fitting planes and interpolating, make sure that the chosen
    %normal makes sense. Particularly, if the chosen normal is along the
    %direction of largest variation, then the corresponding eigenvalue
    %should generally increase, reflecting anaphase behavior. This helps
    %avoid strange prophase/early prometaphase cases which can be slightly
    %elongated and so a plane gets fitted perpendicular to the long
    %direction (in this case, the eigenvalue does not increase over time).
    
    %get the eigenvalues corresponding to the chosen normal and to the
    %other direction(s)
    eigenValOrdered = NaN*(ones(probDim,nTimePoints));
    for t = 1 : length(goodFrames)
        eigenValOrdered(:,goodFrames(t)) = eigenValues(goodFrames(t),...
            eigenVecAssign(:,t));
    end
    eigenValOrdered = [eigenValOrdered(normalIndx,:); ...
        eigenValOrdered(setxor((1:probDim),normalIndx),:)];
    
    %find frames where the first eigenvalue is largest
    frames2consider = [];
    for t = goodFrames
        if eigenValOrdered(1,t) == max(eigenValOrdered(:,t))
            frames2consider = [frames2consider t];
        end
    end

    if length(frames2consider) > 1

        %calculate the eigenvalue ratio
        eigenValRatio = eigenValOrdered(1,frames2consider) ...
            ./ mean(eigenValOrdered(2:probDim,frames2consider));

        %among the frames to consider, find the frames where the eigenratio is
        %larger than the minimum
        framesLarger = find( eigenValRatio >= minEigenValRatio2 );
        
        %do not fit planes in frames where the eigenValRatio is smaller
        %than the minimum
        goodFrames = setdiff(goodFrames,frames2consider(eigenValRatio < ...
            minEigenValRatio2));

        %if 2 or more frames satisfy this criterion, fit a straight line in the
        %eigenvalue ratio and proceed with plane fitting in these frames
        %only if the line has a significant positive slope
        if length(framesLarger) >= 2

            %fit straight line
            [lineFitParam,S] = polyfit(framesLarger,eigenValRatio(framesLarger),1);

            %get std of slope
            varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
            slopeStd = sqrt(varCovMat(1,1));
            
            %test significance of slope and decide whether to proceed with
            %plane fitting
            testStat = lineFitParam(1)/slopeStd;
            pValue = 1 - tcdf(testStat,S.df);
            if pValue > 0.0001
                goodFrames = setdiff(goodFrames,frames2consider(framesLarger));
            end

        end %(if length(framesLarger) >= 2)

    end %(if length(frames2consider) > 1)

else

    %go over frames without a plane and make all kinetochores inliers
    for t = 1 : nTimePoints
        planeFit(t).inlierIdx = (1:nSpots(t));
        planeFit(t).unalignedIdx = [];
        planeFit(t).laggingIdx = [];
    end
    
end %(if nConsecFrames >= minConsecFrames)

%Put back the third dimension.
%For now, this is done in a simple way, where the metaphase plate
%is assumed to be absolutely perpendicular to the x,y-plane.
%This can be improved in the future by rotating the normal to
%make an angle with the x,y-plane such that the positional scatter
%along the normal direction is minimized. - KJ
if use2D
    for t = 1 : nTimePoints
        eigenVecTmp = planeFit(goodFrames(t)).eigenVectors;
        planeFit(goodFrames(t)).eigenVectors = [[eigenVecTmp; zeros(1,2)] [0 0 1]'];
        planeFit(goodFrames(t)).eigenValues = [planeFit(goodFrames(t)).eigenValues NaN];
    end
    eigenVectorsOld = eigenVectors;
    eigenVectors = zeros(3,3,nTimePoints);
    eigenVectors(1:2,1:2,:) = eigenVectorsOld;
    eigenVectors(3,3,:) = 1;
end

if nConsecFrames >= minConsecFrames && ~isempty(goodFrames)

    for t = 1:length(goodFrames)

        % define plane vectors etc.
        goodNormals(:,t) = eigenVectors(:,eigenVecAssign(normalIndx,t),goodFrames(t));
        e_plane = calcPlaneVectors(goodNormals(:,t));
        planeFit(goodFrames(t)).plane = [goodNormals(:,t)',meanCoordFull(goodFrames(t),:)*goodNormals(:,t)];
        planeFit(goodFrames(t)).planeVectors = e_plane;

        % assignment of metaphase or anaphase; distinction of late prometaphase
        % to metaphase will be done during the search for unaligned
        % kinetochores
        eigenRatio = eigenValues(goodFrames(t),eigenVecAssign(normalIndx,t))/...
            mean(eigenValues(goodFrames(t),setdiff(1:probDim,eigenVecAssign(normalIndx,t))));
        if(eigenRatio > 1)
            planeFit(goodFrames(t)).phase = 'a';
        else
            planeFit(goodFrames(t)).phase = 'm';
        end
        
    end

    % find frames whose normal needs to be interpolated; no extrapolation is
    % being done
    gapFrames = setdiff(goodFrames(1):goodFrames(end),goodFrames);

    if ~isempty(gapFrames)
        
        %         % B-spline interpolation of the normals in gap frames
        %         % The interpolation is forced to use derivative 0 at the boundary frames
        %         gapNormals = spline(goodFrames,[[0;0;0] goodNormals [0;0;0]], gapFrames);

        % an attempt to check out linear interpolation - seems to work
        % better than spline
        gapNormals = interp1q(goodFrames',goodNormals',gapFrames');
        gapNormals = gapNormals';
        
        % normalization of interpolated vectors
        gapNormalsNorm = sqrt(sum(gapNormals.^2));
        for i = 1:size(gapNormals,2)
            gapNormals(:,i) = gapNormals(:,i)/gapNormalsNorm(i);
        end

        for t = 1:length(gapFrames)

            % define the interpolated plane vectors etc.
            e_plane = calcPlaneVectors(gapNormals(:,t));
            planeFit(gapFrames(t)).plane = [gapNormals(:,t)',meanCoord(gapFrames(t),:)*gapNormals(:,t)];
            planeFit(gapFrames(t)).planeVectors = e_plane;

            % assign the mitotic phase to what the next good frame is classified
            % as, i.e. in a sequence 'm' 'm' 'e' 'm', 'e' will be replaced by 'm';
            % in a sequence 'm' 'm' 'e' 'e' 'e' 'a' 'a', all 'e's will be replaced
            % by 'a's;
            % dependent on the noise level there might be inconsistencies in the
            % classification. Those will be fetched later in a global consistency
            % check of the mitotic phase classification
            nextGoodFrames = find(goodFrames > gapFrames(t));
            planeFit(gapFrames(t)).phase = planeFit(goodFrames(nextGoodFrames(1))).phase;
            
        end
        
    end %(if ~isempty(gapFrames))

    % get all classified frames (goodFrames and gapFrames)
    framesWiPlane = (goodFrames(1):goodFrames(end));
    framesNoPlane = setxor((1:nTimePoints),framesWiPlane);

    % refinement of phase classification based on scatter along normal to
    % plane and time evolution
    
    % get distance from plane and in-plane coordinates by transformation
    % with inverse of in-plane vectors
    dist2planeStd = NaN(nTimePoints,1);
    for t = framesWiPlane(1):framesWiPlane(end)

        % planeCoord: [d,xplane,yplane]
        planeFit(t).planeCoord = ...
            (inv(planeFit(t).planeVectors)*...
            (initCoord(t).allCoord(:,1:3)-...
            repmat(planeFit(t).planeOrigin,nSpots(t),1))')';
        
        % calculate std of distances from plane
        dist2planeStd(t) = std(planeFit(t).planeCoord(planeFit(t).inlierIdx,1));

    end
    
    %find last frame not labeled 'e' and first frame labeled 'a'
    framePhase = vertcat(planeFit.phase);
    lastFrameNotE = find(framePhase~='e',1,'last');
    firstFrameA = find(framePhase=='a',1,'first');

    %if there is a frame labeled 'a', go back in time and, as long as the
    %kinetochore scatter decreases, classify preceding frames as 'a'
    if ~isempty(firstFrameA)

        t = firstFrameA - 1;
        while t > 0 && (all(dist2planeStd(t) > dist2planeStd(max(1,t-5):max(1,t-1))))
            planeFit(t).phase = 'a';
            t = t - 1;
        end

        %if there are no 'a' frames but there are some empty frames toward
        %the end of the movie, label those as 'a' and again go back in time
        %looking for the start of anaphase
    elseif ~isempty(lastFrameNotE) && lastFrameNotE < nTimePoints
        for t = lastFrameNotE+1 : nTimePoints
            planeFit(t).phase = 'a';
        end
        t = lastFrameNotE;
        while t > 0 && (all(dist2planeStd(t) > dist2planeStd(max(1,t-5):max(1,t-1))))
            planeFit(t).phase = 'a';
            t = t - 1;
        end

    end
        
    % identification of unaligned and lagging kinetochores and further
    % refinement of phase classification

    % identify outliers in each frame
    for t = framesWiPlane(1):framesWiPlane(end)

        %extract distance from plane along the normal
        d = planeFit(t).planeCoord(:,1);

        %if this is a 'm' frame
        if strcmp(planeFit(t).phase,'m') || strcmp(planeFit(t).phase,'p')

            %detect outliers based on d
            [outlierIdx,inlierIdx] = detectOutliers(d,2.5);

            %if there are outliers ...
            if ~isempty(outlierIdx)

                %change label to 'l' ( = late prometaphase)
                planeFit(t).phase = 'p';

                %put indices in their place
                planeFit(t).unalignedIdx = outlierIdx';
                planeFit(t).inlierIdx = inlierIdx';

            else

                %change label to 'm' ( = metaphase)
                planeFit(t).phase = 'm';

                %indicate that there are no outliers
                planeFit(t).unalignedIdx = [];
                planeFit(t).inlierIdx = inlierIdx';

            end

        elseif strcmp(planeFit(t).phase,'a')

            %put distances into two groups: +ve and -ve distances from the
            %plane
            indxNeg = find(d <= 0);
            distNeg = d(indxNeg);
            indxPos = find(d > 0);
            distPos = d(indxPos);

            %find outliers in each group
            [outlierIdxNeg,inlierIdxNeg] = detectOutliers(distNeg,3);
            [outlierIdxPos,inlierIdxPos] = detectOutliers(distPos,3);

            %add outliers ahead of the majority to the list of unaligned
            %kinetochores
            %add outliers lagging behind the majority to the list of
            %lagging kinetochores
            unalignedIdx = [];
            laggingIdx = [];
            if ~isempty(outlierIdxNeg)
                distOutlier = distNeg(outlierIdxNeg);
                meanNeg = mean(distNeg(inlierIdxNeg));
                unalignedIdx = (indxNeg(outlierIdxNeg(...
                    distOutlier <= meanNeg)))';
                laggingIdx = (indxNeg(outlierIdxNeg(...
                    distOutlier > meanNeg)))';
            end
            if ~isempty(outlierIdxPos)
                distOutlier = distPos(outlierIdxPos);
                meanPos = mean(distPos(inlierIdxPos));
                unalignedIdx = [unalignedIdx ...
                    (indxPos(outlierIdxPos(distOutlier >= meanPos)))'];
                laggingIdx = [laggingIdx ...
                    (indxPos(outlierIdxPos(distOutlier < meanPos)))'];
            end

            %save indices in planeFit
            if isempty(unalignedIdx)
                planeFit(t).unalignedIdx = [];
            else
                planeFit(t).unalignedIdx = unalignedIdx;
            end
            if isempty(laggingIdx)
                planeFit(t).laggingIdx = [];
            else
                planeFit(t).laggingIdx = laggingIdx;
            end
            planeFit(t).inlierIdx = [indxPos(inlierIdxPos)' indxNeg(inlierIdxNeg)'];

        end %(if strcmp(planeFit(t).phase,'m') ...)

    end %(for t = framesWiPlane(1):framesWiPlane(end))
    
    %go over frames without a plane and make all kinetochores inliers
    for t = framesNoPlane
        planeFit(t).inlierIdx = (1:nSpots(t));
        planeFit(t).unalignedIdx = [];
        planeFit(t).laggingIdx = [];        
    end

else
    
    %go over frames without a plane and make all kinetochores inliers
    for t = 1 : nTimePoints
        planeFit(t).inlierIdx = (1:nSpots(t));
        planeFit(t).unalignedIdx = [];
        planeFit(t).laggingIdx = [];        
    end    

end %(if nConsecFrames >= minConsecFrames && ~isempty(goodFrames))

%make frame classification (e,p,m,a) smooth, i.e. there should be consecutive
%sequences of e's, then p's, then m's and then a's, with no intermingling

%get frame phases
framePhase = vertcat(planeFit.phase);

%convert letters to numbers
framePhaseNum = zeros(size(framePhase));
framePhaseNum(framePhase=='p') = 1;
framePhaseNum(framePhase=='m') = 2;
framePhaseNum(framePhase=='a') = 3;

%find last frame that has a nonzero phase to start reassignment from 
lastFrame = find(framePhaseNum~=0,1,'last');

%go backwards over all frames and assign correct phase
for t = lastFrame-1 : -1 : 1
    framePhaseNum(t) = min(framePhaseNum(t:lastFrame));
end

%convert back to letters
framePhase(framePhaseNum==0) = 'e';
framePhase(framePhaseNum==1) = 'p';
framePhase(framePhaseNum==2) = 'm';
framePhase(framePhaseNum==3) = 'a';

%put phases back in structure
for t = 1 : nTimePoints
    planeFit(t).phase = framePhase(t);
end

%% align frames wherever possible to get rid of overall rotation

%for frames with a plane, use plane origin as the frame origin
%for frames without a plane, use the center of mass as the frame origin
frameOrigin = vertcat(planeFit.planeOrigin);

%shift the coordinates in each frame such that frameOrigin 
%in each frame is the origin
tmpCoord = repmat(struct('allCoord',[]),nTimePoints,1);
for iTime = 1 : nTimePoints
    tmpCoord(iTime).allCoord = initCoord(iTime).allCoord;
    tmpCoord(iTime).allCoord(:,1:3) = tmpCoord(iTime).allCoord(:,1:3) - ...
        repmat(frameOrigin(iTime,:),nSpots(iTime),1);
end

%if there are frames to rotate ...
if length(framesWiPlane) > 1

    %get first frame to rotate
    firstFrameRotate = framesWiPlane(1);

    %find frames without plane that are before the first frame with plane
    framesBefore1 = framesNoPlane(framesNoPlane < firstFrameRotate);
    framesAfter1 = setxor(framesNoPlane,framesBefore1); % the rest of the frames

    %get the coordinate system of each frame with a plane
    coordSystem = zeros(3,3,nTimePoints);
    coordSystem(:,:,framesWiPlane) = cat(3,planeFit.planeVectors);

    %assign the coordinate system of each frame without a plane
    %if they are before the first frame with a plane, take the coordinate
    %system of the plane after
    %if they are after the first frame with a plane, take the coordinate
    %system of the plane before
    for iTime = framesBefore1(end:-1:1)
        coordSystem(:,:,iTime) = coordSystem(:,:,iTime+1);
    end
    for iTime = framesAfter1
        coordSystem(:,:,iTime) = coordSystem(:,:,iTime-1);
    end
    
    %rotate the coordinates in all frames
    %propagate errors to the new coordinates
    for iTime = 1 : nTimePoints
        rotationMat = inv(coordSystem(:,:,iTime)); %rotation matrix
        errorPropMat = rotationMat.^2; %error propagation matrix
        tmpCoord(iTime).allCoord(:,1:3) = (rotationMat*(tmpCoord(iTime).allCoord(:,1:3))')';
        tmpCoord(iTime).allCoord(:,4:6) = sqrt((errorPropMat*((tmpCoord(iTime).allCoord(:,4:6)).^2)')');
    end
    
end %(if length(framesWiPlane) > 1)

%store rotated coordinates in planeFit structure
for iTime = 1 : nTimePoints
    planeFit(iTime).rotatedCoord = tmpCoord(iTime).allCoord;
end

%% output

% assign output
dataStruct.planeFit = planeFit;

% plot everything if verbose
if verbose > 0
    % makiFitPlanePlot(dataStruct)
    makiPlotRotatingPlanes(dataStruct);
end

% turn warnings back on
warning(warningState);



%% SUBFUNCTIONS

function [eigenVectors, eigenValues, meanCoord] = eigenCalc(coordinates)

%get problem dimensionality
probDim = size(coordinates,2);

% find eigenVectors, eigenValues, meanCoordinates for plane fitting
coordCov = cov(coordinates);

meanCoord = mean(coordinates,1);

% eigenVectors/eigenValues
[eigenVectors,eVals] = eig(coordCov);
eigenValues = diag(eVals)';

if probDim==3

    % compare eigenValues
    diffs = pdist(eigenValues');

    % indices to understand lowest
    [u,v] = find(tril(ones(probDim),-1));
    [dummy,idx] = min(diffs);

    % find two close, one far
    closestIdx = [u(idx),v(idx)];
    farIdx = setdiff(1:probDim,closestIdx);
    
else
   
    %The algorithm above does not work in the case of 2D. Take simple
    %alternative, which makes use of the fact that the first eigenvalue is
    %always smallest. Assumes metaphase/late prometaphase configuration.    
    farIdx = 1;
    closestIdx = 2;
    
end


% sort eigenVectors, eigenValues. X-component of first eigenvector should
% always be positive
eigenVectors = eigenVectors(:,[farIdx,closestIdx]).*sign(eigenVectors(1,farIdx));
eigenValues = eigenValues([farIdx,closestIdx]);


function [aList, evec, rotList] = assignEigenVecs(evec) 
% assign the indices of the eignevectors such that the tripod between
% consecutive frames undergoes minimal overall rotation
% The output of the function is a probDim x nTimePoints matrix where columns list
% eigenvector index for a specific frame, while the next column lists the
% indices of the associated eigenvector in the next frame. For instance, in
% 3D,
%           1 1 1 1 2 1 1 ...
% aList = [ 2 2 3 2 1 3 3 ...
%           3 3 2 3 3 2 2 ...
%
% indicates that the eigenvectors in order of input have the best match
% between frame 1 and 2; but between frame 2 and 3 eigenvectors 2 and 3
% have been swapped, etc.
%
% The function also updates the eigenvectors to remove near-180 degrees
% jumps
%
% rotList contains the cos(Phi_t->t+1) for each the eigenvector assignments

[probDim,probDim,nTimePoints] = size(evec);
aList = zeros(probDim,nTimePoints);
aList(:,1) = (1:probDim)'; 

for i = 2:nTimePoints
    % calculate the cost matrix eignevector assignment
    costMat = transpose(evec(:,:,i-1))*evec(:,:,i);
    % the cost for a rotation is 1 - abs(cos(evec_i(t)*evec_j(t+1))
    costMat = 1 - abs(costMat);
    links = lap(costMat);
    aList(:,i)=links(aList(:,i-1));
    rotList(:,i-1)= diag(1 - costMat(aList(:,i-1),aList(:,i)));
    
    % check for negative dot products -- these indicate 180 degree jumps
    for j = 1:probDim
        if(transpose(evec(:,aList(j,i-1),i-1))*evec(:,aList(j,i),i) < 0)
            evec(:,aList(j,i),i)=-1*evec(:,aList(j,i),i);
        end
    end
end

function e_plane = calcPlaneVectors(normal)
% calculates the plane vectors from a normal based on the criterion that
% the first plane vector is the normal, the second is the vector
% perpependicular to the normal and parallel to the XY-plane and the third
% vector is perpendicular to both

e_plane = zeros(3);
e_plane(:,1) = normal;
e_plane(:,2) = [-normal(2),normal(1),0]./sqrt(sum(normal(1:2).^2));
e_plane(:,3) = cross(e_plane(:,1),e_plane(:,2));



%% old stuff

    % % do only for good frames
    % if ~isempty(metaphaseFrames)
    %     for t = metaphaseFrames'
    %         done = false;
    %         % initially: assume no bad spots. Allow for 10 iterations
    %         badSpotIdxLOld = false(nSpots(t),10);
    %         ct = 1;
    %
    %         while ~done
    %
    %             % get distance from plane, in-plane coordinates by transformation
    %             % with inverse of in-plane vectors
    %             normal = eigenVectors(:,1,t);
    %             e_plane = zeros(3);
    %             e_plane(:,1) = normal;
    %             e_plane(:,2) = [-normal(2),normal(1),0]./sqrt(sum(normal(1:2).^2));
    %             e_plane(:,3) = cross(e_plane(:,1),e_plane(:,2));
    %             % planeCoord: [d,xplane,yplane]
    %             planeFit(t).planeCoord = ...
    %                 (inv(e_plane)*...
    %                 (initCoord(t).allCoord(:,1:3)-repmat(meanCoord(t,:),nSpots(t),1))')';
    %
    %             % find outliers
    %             [dummy, dummy, goodSpotIdx] = robustMean(planeFit(t).planeCoord(:,1));
    %             badSpotIdxL = true(initCoord(t).nSpots,1);
    %             badSpotIdxL(goodSpotIdx) = false;
    %
    %             %         % DEBUG plot plane in matlab
    %             %         pc=planeFit(t).planeCoord;
    %             %         pos = pc(:,1)>0;
    %             %         neg = pc(:,1)<0;
    %             %         figure,plot3(pc(pos,1),pc(pos,2),pc(pos,3),'.k',...
    %             %             pc(neg,1),pc(neg,2),pc(neg,3),'or',...
    %             %             pc(badSpotIdxL,1),pc(badSpotIdxL,2),pc(badSpotIdxL,3),'+b')
    %
    %             % check whether there was any change
    %             if any(all(repmat(badSpotIdxL,1,10) == badSpotIdxLOld,1)) || ct == 10
    %                 % done. Fill information into planeFit-structure
    %
    %                 planeFit(t).plane = [normal',meanCoord(t,:)*normal];
    %                 planeFit(t).planeVectors = e_plane;
    %                 planeFit(t).eigenVectors = eigenVectors(:,:,t);
    %                 planeFit(t).eigenValues = eigenValues(t,:);
    %                 planeFit(t).planeOrigin = meanCoord(t,:);
    %
    %                 % lagging chromosomes are outliers (until we can identify pairs
    %                 planeFit(t).laggingIdx = find(badSpotIdxL);
    %                 planeFit(t).inlierIdx = goodSpotIdx;
    %
    %                 % distribution parameters (do for all unit directions -
    %                 % the second vector is also interesting, as it lies in the xy
    %                 % plane, in which the metaphase plate should not be cut off
    %                 % distribution parameters (rows):
    %                 % var
    %                 % skew
    %                 % kurtosis
    %                 % p for normal distribution (lilliefors test)
    %                 % correct all parameters for bias
    %                 planeFit(t).distParms = zeros(4,3);
    %                 planeFit(t).distParms(1,:) = var(planeFit(t).planeCoord(goodSpotIdx,:));
    %                 planeFit(t).distParms(2,:) = skewness(planeFit(t).planeCoord(goodSpotIdx,:),0);
    %                 planeFit(t).distParms(3,:) = kurtosis(planeFit(t).planeCoord(goodSpotIdx,:),0);
    %                 [dummy,planeFit(t).distParms(4,1)] = ...
    %                     lillietest(planeFit(t).planeCoord(goodSpotIdx,1));
    %                 [dummy,planeFit(t).distParms(4,2)] = ...
    %                     lillietest(planeFit(t).planeCoord(goodSpotIdx,2));
    %                 [dummy,planeFit(t).distParms(4,3)] = ...
    %                     lillietest(planeFit(t).planeCoord(goodSpotIdx,3));
    %
    %                 % plot plane in matlab
    %                 if verbose > 1
    %                     [ygrid,zgrid] = meshgrid(...
    %                         linspace(min(planeFit(t).planeCoord(:,2)),...
    %                         max(planeFit(t).planeCoord(:,2)),5), ...
    %                         linspace(min(planeFit(t).planeCoord(:,3)),...
    %                         max(planeFit(t).planeCoord(:,3)),5));
    %                     xgrid = zeros(5,5);
    %                     pc=planeFit(t).planeCoord;
    %                     pos = pc(:,1)>0;
    %                     neg = pc(:,1)<0;
    %                     figure('Name',...
    %                         sprintf('Metaphase plate frame %i for %s',...
    %                         t,dataStruct.projectName))
    %                     plot3(pc(pos,1),pc(pos,2),pc(pos,3),'.k',...
    %                         pc(neg,1),pc(neg,2),pc(neg,3),'or',...
    %                         pc(badSpotIdxL,1),pc(badSpotIdxL,2),pc(badSpotIdxL,3),'+b')
    %                     hold on
    %                     mesh(xgrid,ygrid,zgrid,'EdgeColor',[0 0 1],'FaceAlpha',0);
    %                     grid on
    %                 end
    %
    %
    %                 done = true;
    %
    %             else
    %                 % re-"fit" the plane. Update eigenVectors etc.
    %                 [eigenVectors(:,:,t), eigenValues(t,:), meanCoord(t,:)] = ...
    %                     eigenCalc(initCoord(t).allCoord(goodSpotIdx,1:3));
    %                 % count fit
    %                 ct = ct + 1;
    %                 % remember current bad spots
    %                 badSpotIdxLOld(:,ct) = badSpotIdxL;
    %
    %             end
    %
    %         end
    %     end % loop good frames
    % end

    % % loop to get the "between frames" - stuff
    % if length(metaphaseFrames) > 1
    %     for t=[metaphaseFrames(1:end-1),metaphaseFrames(2:end)]'
    %         % p-value of distribution comparison
    %         [dummy,planeFit(t(1)).deltaP] = kstest2(...
    %             planeFit(t(1)).planeCoord(planeFit(t(1)).inlierIdx,1),...
    %             planeFit(t(2)).planeCoord(planeFit(t(2)).inlierIdx,1));
    %
    %         % change in plane orientation
    %         planeFit(t(1)).deltaAngle = acos(dot(planeFit(t(1)).planeVectors(:,1),...
    %             planeFit(t(2)).planeVectors(:,1))) *180/pi;
    %     end
    % end

    
            %         %construct a hierarchical clustering tree with complete linkage
        %         cTree = linkage(pdist(d),'average');
        %         dendrogram(cTree,0);
        %         input(sprintf('frame %3d: enter to continue',t));
        %
        %         %find groups of kinetochores based on this clustering
        %         kinetLabels = cluster(cTree,'cutoff',1,'criterion','inconsistent');
        %         numGroups = max(kinetLabels);
        %         kinetGroups = repmat(struct('indx',[],'distance',[]),numGroups,1);
        %         distMean = zeros(numGroups,1);
        %         for iGroup = 1 : numGroups
        %             kinetGroups(iGroup).indx = find(kinetLabels==iGroup);
        %             kinetGroups(iGroup).distance = d(kinetGroups(iGroup).indx);
        %             distMean(iGroup) = mean(kinetGroups(iGroup).distance);
        %         end
        %         distMean = abs(distMean);
        %
        %         %identify inlier and outlier indices
        %         %inlier group is the one with smallest mean distance
        %         %all other groups are outliers
        %         inlierGroup = find(distMean==min(distMean));
        %         inlierIdx = kinetGroups(inlierGroup).indx;
        %         outlierGroup = setxor((1:numGroups),inlierGroup);
        %         outlierIdx = vertcat(kinetGroups(outlierGroup).indx);

