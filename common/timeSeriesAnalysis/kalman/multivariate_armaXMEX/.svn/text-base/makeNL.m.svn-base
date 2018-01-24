matOpsDir = ' ./Utilities/Matrix_Operations/';
utlDir = ' ./Utilities/Carma_Utilities/';
minDir = ' ./Minimizer/';
carmaDir = ' ./';

%flags = ' -g -output carmaFitModel '; %For debugging
flags = ' -output carmaFitModel '; %for speed (with optimizations)

callString = 'mex';

matOpsArray = {'matrixAdd.c',
               'matrixMultiply.c',
               'matrixMultiplyConst.c',
               'matrixSubtract.c',
               'matrixTranspose.c',
               'leftdiv.c',
               'matOps.h'};
           
utilsArray = {'inverseLevinsonDurbinExpoAR.c',
              'inverseLevinsonDurbinExpoMA.c',
              'levinsonDurbinExpoAR.c',
              'levinsonDurbinExpoMA.c',
              'paramsFromVector.c',
              'vectorFromParams.c',
              'carmaUtils.h'};

carmaArray = {'writeLogFile.c',
              'errHandle.c',
              'carmaFitModel.c',
              'carmaCalcKalmanInnov.c',
              'covKalmanInit.c',
              'carmaNegLnLikelihood.c',
              'prob.c',
              'prob.h',              
              'carmaFitModel.h'};
          
minimArray = {'amoeba.c',
              'amoeba.h',
              'nrutil.c',              
              'nrutil.h',
              'createSimplex.c'};
          
matOpsString = '';

for j=1:length(matOpsArray)
    
  matOpsString = strcat(matOpsString, matOpsDir);
  matOpsString = strcat(matOpsString, cast(matOpsArray(j),'char'));

end


utilsString = '';

for j=1:length(utilsArray)
    
  utilsString = strcat(utilsString, utlDir);
  utilsString = strcat(utilsString, cast(utilsArray(j),'char'));

end

carmaString = '';

for j=1:length(carmaArray)
  
  carmaString = strcat(carmaString, carmaDir);
  carmaString = strcat(carmaString, cast(carmaArray(j),'char'));

end

minimString = '';

for j=1:length(minimArray)
  
  minimString = strcat(minimString, minDir);
  minimString = strcat(minimString, cast(minimArray(j),'char'));

end




callString = strcat(callString,flags);
callString = strcat(callString,carmaString);
callString = strcat(callString,utilsString);
callString = strcat(callString,matOpsString);
callString = strcat(callString,minimString);

eval(callString);

clear j callString flags carmaString utilsString matOpsString minimString
clear matOpsDir utlDir carmaDir minDir 
clear matOpsArray utilsArray carmaArray minimArray
