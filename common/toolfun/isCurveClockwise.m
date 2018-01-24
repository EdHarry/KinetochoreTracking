function isClockwise = isCurveClockwise(curveIn)

%Determines wether a (closed) input curve runs clockwise or counter
%clockwise

%WARNING: when used with really wierd shaped curves, the spline tolerance
%may need to be increased to get accurate results

isClockwise = [];

if ndims(curveIn) ~= 2
    error('input curve of wrong dimension!');
    return
end

if length(curveIn) < 3
    error('Need at least three points to make a closed curve!')
    return
end

%Transpose to row vector if necessary
if size(curveIn,1) > size(curveIn,2)
    curveIn = curveIn';
end

%use the polygon area formula to check curve orientation
v1 = curveIn(1,1:end-1) .* curveIn(2,2:end);
v2 = curveIn(1,2:end) .* curveIn(2,1:end-1);

isClockwise = (sum(v1 - v2) / 2) < 0 ;