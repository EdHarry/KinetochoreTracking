function fitRes = fit1Line2Lines(xData,yData,alpha)

%OUTPUT fitRes: Structure with fields:
%           .lines1or2: Better fit (1 or 2 lines).
%           .line1: Parameters of 1 line fit: 
%               .param: slope and y-intercept.
%               .residStd: residuals std.
%               .varCovMat: variance-covariance matrix of parameters.
%           .line2: Parameters of 2 lines fit: 
%               .param: 1st line slope, 1st line y-intercept, 2nd line
%                   slope, 2nd line y-intercept and bounrary dividing the
%                   two regimes.
%               .residStd: residuals std.
%               .varCovMat: variance-covariance matrix of parameters.

options = optimset('Display','off');

%fit data with one line
initialGuess = [1 0]';
[param1Line,resnorm,resid1Line,dummy,dummy,dummy,jac1Line] = ...
    lsqcurvefit(@fit1line,initialGuess,xData,yData,[],[],options);
numDegFree1Line = length(xData) - 2;
varCov1Line = full(inv(jac1Line' * jac1Line) * (resid1Line' * ...
    resid1Line) / numDegFree1Line);

%fit data with two lines
initialGuess = [param1Line; param1Line(1); 0.5];
[param2Lines,resnorm,resid2Lines,dummy,dummy,dummy,jac2Lines] = ...
    lsqcurvefit(@fit2lines,initialGuess,xData,yData,[],[],options);
numDegFree2Lines = length(xData) - 4;
varCov2LinesTmp = full(inv(jac2Lines' * jac2Lines) * (resid2Lines' * ...
    resid2Lines) / numDegFree2Lines);
varCov2Lines = NaN(5);
varCov2Lines(1:3,1:3) = varCov2LinesTmp(1:3,1:3);
varCov2Lines(5,1:3) = varCov2LinesTmp(4,1:3);
varCov2Lines(5,5) = varCov2LinesTmp(4,4);
varCov2Lines(1:3,5) = varCov2LinesTmp(1:3,4);

%compare residuals to decide which is the better fit
%get test statistic, which is F-distributed
testStat = (sum(resid2Lines.^2)/numDegFree2Lines)/...
    (sum(resid1Line.^2)/numDegFree1Line);

%get p-value of test statistic
pValue = fcdf(testStat,numDegFree2Lines,numDegFree1Line);

%compare p-value to alpha-value (1-sided F-test)
%if pValue < alpha-value, 2 lines; otherwise, 1 line
numLines = (pValue < alpha) + 1;

%output
fitRes.lines1or2 = numLines;
fitRes.line1.param = param1Line;
fitRes.line1.residStd = nanstd(resid1Line);
fitRes.line1.varCovMat = varCov1Line;
fitRes.line2.param = [param2Lines(1:3); (param2Lines(1) - param2Lines(3)) * ...
    param2Lines(4) + param2Lines(2); param2Lines(4)];
fitRes.line2.residStd = nanstd(resid2Lines);
fitRes.line2.varCovMat = varCov2Lines;

function yData1Line = fit1line(param,xData)

yData1Line = param(1) * xData + param(2);
yData1Line = yData1Line(~isnan(yData1Line));

function yData2Lines = fit2lines(param,xData)

yData2Lines = NaN(size(xData));
yData2Lines(xData<param(4)) = param(1) * xData(xData<param(4)) + param(2);
yData2Lines(xData>=param(4)) = param(3) * xData(xData>=param(4)) + ...
    (param(1) - param(3)) * param(4) + param(2);
yData2Lines = yData2Lines(~isnan(yData2Lines));
