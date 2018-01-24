function [ timeWindowFinal,minRadiusFinal,maxRadiusFinal ] = trackingOptimisation( job )
% EHarry Jan 2012

x0=[20,0.1,0.2,0.1,0.5,1.5,0.5]./1e6;
%x0 = [job(1).dataStruct.dataProperties.tracksParam.costMatrices(1).parameters.nnWindow-1,job(1).dataStruct.dataProperties.tracksParam.costMatrices(1).parameters.minSearchRadius,job(1).dataStruct.dataProperties.tracksParam.costMatrices(1).parameters.maxSearchRadius]./1e10;
% lb = [0,0,0,0,0.1,0.1,0.1];
% ub = [30,3,3,3,10,10,10];

%x0 = [];

op = optimset('Display','off','MaxIter',10000,'MaxFunEvals',10000,'TolFun',1e-10,'TolX',1e-10);

parametersFinal = lsqnonlin(@(parameters)optimisation(job,parameters),x0,[],[],op);

parametersFinal = parametersFinal .* 1e6;

timeWindowFinal = round(abs(parametersFinal(1)));
minRadiusFinal = abs(parametersFinal(2:4));
maxRadiusFinal = abs(parametersFinal(5:7));

%% SUBFUNCTIONS

    function dataStruct = generateTrackingParam(dataStruct,timeWindowTmp,minRadiusTmp,maxRadiusTmp)
        
        
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


    function cost = optimisation(job,parameters)
        
        parameters = parameters .* 1e6;
        
        timeWindowTmp = round(abs(parameters(1)));
        minRadiusTmp = abs(parameters(2:4));
        maxRadiusTmp = abs(parameters(5:7));
        
        job = makiMakeJob_poleOME('MCAINSH',[1 2 3 5 -6 -7 -8 -9 -10 -11 -12 -13],job,0);
        %  job = makiMakeJobPlatformIndependent(job,'MCAINSH');
        
        for i = 1:length(job)
            job(i).dataStruct = generateTrackingParam(job(i).dataStruct,timeWindowTmp,minRadiusTmp,maxRadiusTmp);
        end
        
        display(['Current Parameters: ',num2str(parameters)]);
        
        %makiMovieAnalysis_poleOME('MCAINSH',job);
        
        master_track_poleOME(job,[],0,0);
        
       % job = makiMakeJobPlatformIndependent(job,'MCAINSH');
        job = makiMakeJob_poleOME('MCAINSH',[],job,0);
        
        cost = fitTracking( job );
        
        display(['Current Cost: ',num2str(cost)]);
    end

end

