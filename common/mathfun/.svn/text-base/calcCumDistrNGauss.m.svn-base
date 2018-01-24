function cumDistrNGauss = calcCumDistrNGauss(param,abscissa,variableMean,...
    variableStd)
%CALCCUMDISTRNGAUSS calculates the cumulative distribution of N Gaussians
%
%SYNOPSIS cumDistrNGauss = calcCumDistrNGauss(param,abscissa,variableMean,...
%    variableStd)
%
%INPUT  param         : Vector of parameters indicating the means,
%                       variances and amplitudes of the N Gaussians.
%                       -If variableMean=1 & variableStd=1, param has 3N 
%                        entries: N means, N standard deviations and N amplitudes.
%                       -If variableMean=1 & variableStd=0, param has 2N+1
%                        entries: N means, 1 standard deviation and N amplitudes.
%                       -If variableMean=1 & variableStd=2, param has 2N+1
%                        entries: N means, 1 standard deviation and N amplitudes.
%                       -If variableMean=0 & variableStd=1, param has 2N+1
%                        entries: 1 mean, N standard deviations and N amplitudes.
%                       -If variableMean=0 & variableStd=0, param has N+2
%                        entries: 1 mean, 1 standard deviation and N amplitudes.
%                       See below the definitions of variableMean and
%                       variableStd.
%       abscissa      : Abscissa values at which the cumulative
%                       distribution is calculated.
%       variableMean  : 0 if assuming the fixed relationship
%                       (mean of nth Gaussian) = n * (mean of 1st Gaussian).
%                       1 if there is no relationship between the means of 
%                       different Gaussians. 
%                       Optional. Default: 1.
%       variableStd   : 0 if assuming that all Gaussians have the same
%                       standard deviation. 1 if there is no relationship
%                       between the standard deviations of different
%                       Gaussians. 2 if assuming the relationship (std of
%                       nth Gaussian) = sqrt(n) * (std of 1st Gaussian).
%                       Optional. Default: 1.
%
%OUTPUT cumDistrNGauss: Values of the resulting cumulative distribution
%                       given the input abscissa values.

%Khuloud Jaqaman, August 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cumDistrNGauss = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 2
    disp('--calcCumDistrNGauss: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

if nargin < 3 || isempty(variableMean)
    variableMean = 1;
end

if nargin < 4 || isempty(variableStd)
    variableStd = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the cumulative distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the means, variances and amplitudes of the Gaussians from the input
%parameter vector
switch variableMean

    case 0 %if mean is not variable

        switch variableStd

            case 0 %if variance is constrained to all variances are equal

                %get number of Gaussians
                numGauss = length(param)-2;

                %get their means, variances and amplitudes
                gaussMean = [1:numGauss]'*param(1);
                gaussStd  = repmat(param(2),numGauss,1);
                gaussAmp  = param(3:end);

            case 1 %if variance is variable

                %get number of Gaussians
                numGauss = floor(length(param)/2);

                %get their means, variances and amplitudes
                gaussMean = [1:numGauss]'*param(1);
                gaussStd  = param(2:numGauss+1);
                gaussAmp  = param(numGauss+2:end);
                
            case 2 %if variance is constrained to variance_n = n*variance_1
                
                %get number of Gaussians
                numGauss = length(param)-2;
                
                %get their means, variances and amplitudes
                gaussMean = (1:numGauss)'*param(1);
                gaussStd  = sqrt(1:numGauss)'*param(2);
                gaussAmp  = param(3:end);

        end %(switch variableStd)

    case 1 %if mean is variable

        switch variableStd

            case 0 %if variance is constrained to all variances are equal

                %get number of Gaussians
                numGauss = floor(length(param)/2);

                %get their means, variances and amplitudes
                gaussMean = param(1:numGauss);
                gaussStd  = repmat(param(numGauss+1),numGauss,1);
                gaussAmp  = param(numGauss+2:end);

            case 1 %if variance is variable

                %get number of Gaussians
                numGauss = length(param)/3;

                %get their means, variances and amplitudes
                gaussMean = param(1:numGauss);
                gaussStd  = param(numGauss+1:2*numGauss);
                gaussAmp  = param(2*numGauss+1:end);

        end %(switch variableStd)

end %(switch variableMean)

%calculate the cumulative distribution
cumDistrNGauss = zeros(size(abscissa));
for i=1:numGauss
    cumDistrNGauss = cumDistrNGauss + gaussAmp(i)*normcdf(abscissa,...
        gaussMean(i),gaussStd(i));
end

%%%%% ~~ the end ~~ %%%%%

