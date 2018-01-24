function dataStruct = generateTrackingParam(dataStruct,timeWindowTmp,minRadiusTmp,maxRadiusTmp)
% EHarry October 2011, taken from a subfunction of makiMakeJob.m by KJaqaman         
                
        %assign whether to use rotated coordinates or not
        rotate = 1;
        
        %assign gap closing parameters
        gapCloseParam.timeWindow = timeWindowTmp + 1;
        gapCloseParam.mergeSplit = 0;
        gapCloseParam.minTrackLen = 1;
        
        %assign cost matrix parameters for linking spots between consecutive
        %frames
        costMatrices(1).funcName = 'makiTrackCostMatLink';
        parameters.linearMotion = 0;
        parameters.minSearchRadius = minRadiusTmp;
        parameters.maxSearchRadius = maxRadiusTmp;
        parameters.brownStdMult = 3.5;
        parameters.useLocalDensity = 1;
        parameters.nnWindow = gapCloseParam.timeWindow;
        costMatrices(1).parameters = parameters;
        clear parameters
        
        %assign cost matrix parameters for closing gaps and (in principle)
        %merging and splitting
        costMatrices(2).funcName = 'makiTrackCostMatCloseGaps';
        parameters.linearMotion = 0;
        parameters.minSearchRadius = minRadiusTmp;
        parameters.maxSearchRadius = maxRadiusTmp;
        parameters.brownStdMult = 3.5*ones(gapCloseParam.timeWindow,1);
        % parameters.timeReachConfB = min(2,gapCloseParam.timeWindow);
        parameters.timeReachConfB = min(1,gapCloseParam.timeWindow);
        parameters.lenForClassify = 10;
        parameters.ampRatioLimit = [0.65 4];
        parameters.useLocalDensity = 1;
        parameters.nnWindow = gapCloseParam.timeWindow;
        parameters.linStdMult = 3.5*ones(gapCloseParam.timeWindow,1);
        parameters.timeReachConfL = 1;
        parameters.maxAngleVV = 45;
        costMatrices(2).parameters = parameters;
        clear parameters
        
        %assign Kalman filter function names
        kalmanFunctions.reserveMem = 'kalmanResMemLM';
        kalmanFunctions.initialize = 'kalmanInitLinearMotion';
        kalmanFunctions.calcGain = 'kalmanGainLinearMotion';
        kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
        
        %save tracking parameters in dataStruct
        tracksParam.rotate = rotate;
        tracksParam.gapCloseParam = gapCloseParam;
        tracksParam.costMatrices = costMatrices;
        tracksParam.kalmanFunctions = kalmanFunctions;
        tracksParam.pole = 1;
        dataStruct.dataProperties.tracksParam = tracksParam;
        
        
    end
