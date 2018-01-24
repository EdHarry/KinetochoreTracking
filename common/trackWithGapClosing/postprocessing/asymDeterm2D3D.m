function [asymParam,asymFlag] = asymDeterm2D3D(positions,alpha)
%ASYMDETERM2D3D estimates the asymmetry in a scatter of positions
%
%SYNPOSIS [asymParam,asymFlag] = asymDeterm2D3D(positions,alpha)
%
%INPUT  positions: n-by-2/3 array of positions (x,y,[z]).
%       alpha    : Alpha-value for determining the threshold.
%                  Can take the values 0.2, 0.1 and 0.05.
%                  Optional. Default: 0.1.
%
%OUTPUT asymParam: Parameter estimating asymmetry of positional scatter.
%       asymFlag : 1 if scatter is found to be asymmetric (given the input
%                  alpha-value), 0 otherwise.
%
%Khuloud Jaqaman, December 2007

%% input

%assign default alpha if not input
if nargin < 2 || isempty(alpha)
    alpha = 0.1;
else
    if ~any(alpha == [0.2 0.1 0.05])
        disp('--asymDeterm2D3D: alpha can take only the values 0.05, 0.1 or 0.2');
        disp('                  assigning default value of 0.1');
        alpha = 0.1;
    end
end

%get problem dimensionality (2D or 3D)
probDim = size(positions,2);

%get number of time points in track
numTimePoints = length(find(~isnan(positions(:,1))));

%determine threshold based on alpha-value and dimensionality
switch probDim
    case 2
        switch alpha
            case 0.2  %80th percentile
                asymThresh = [[NaN NaN 3.5 2 1.5 1.4 1.3 1.2 1.2 1.1 1.1 ...
                    1.1 1.1 1.1 1.1 1.05 1.05 1.05 1.05 1.05]'; ...
                    ones(max(numTimePoints-20,1),1)];
            case 0.1  %90th percentile
                asymThresh = [[NaN NaN 5 2.7 2.1 1.8 1.7 1.6 1.5 1.45 1.45 ...
                    1.4 1.4 1.4 1.4 1.4 1.4 1.35 1.35 1.35]'; ...
                    1.3*ones(max(numTimePoints-20,1),1)];
            case 0.05 %95th percentile
                asymThresh = [[NaN NaN 6.0 3.5 2.5 2.2 1.9 1.9 1.9 1.8 ...
                    1.75 1.74 1.74 1.73 1.7 1.7 1.7 1.7 1.67 1.65]'; ...
                    1.6*ones(max(numTimePoints-20,1),1)];
        end
    case 3
        switch alpha
            case 0.2  %80th percentile
                asymThresh = [[NaN NaN 2.2 1.4 1.2 1 1 0.96 0.94 0.92]'; ...
                    0.9*ones(max(numTimePoints-10,1),1)];
            case 0.1  %90th percentile
                asymThresh = [[NaN NaN 2.9 1.9 1.5 1.4 1.3 1.3 1.2 1.2 ...
                    1.2 1.2 1.2 1.2 1.2]'; 1.1*ones(max(numTimePoints-15,1),1)];
            case 0.05 %95th percentile
                asymThresh = [[NaN NaN 3.7 2.2 1.9 1.6 1.5 1.4 1.4 1.4]'; ...
                    1.35*ones(max(numTimePoints-10,1),1)];
        end
end

if numTimePoints > 2

    %% asymmetry calculation

    %calculate the variance-covariance matrix of positions
    posCov = nancov(positions);

    %get the eigen-values of the variance-covariance matrix
    eigenVal = eig(posCov);

    %calculate some intermediate sums
    doubleSum = 0;
    for i = 1 : probDim - 1
        doubleSum = doubleSum + sum( ( eigenVal(i) - eigenVal(i+1:end) ).^2 );
    end
    singleSum = (probDim - 1) * ( sum(eigenVal) )^2;

    %calculate asymmetry parameter
    asymParam = -log( 1 - doubleSum / singleSum );

    %% comparison with threshold

    asymFlag = asymParam > asymThresh(numTimePoints);
    
else
    
    asymParam = NaN;
    asymFlag = NaN;
    
end

%% ~~~ the end ~~~
