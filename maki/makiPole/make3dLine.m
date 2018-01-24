function line = make3dLine( point1,point2 , noPoints )
%MAKE3DLINE Make a 3D line of points between two points 
% EHarry October 2011

vector = point2 - point1;

vectorBase = vector ./ (noPoints-1);

% line starts at point1 and ends at point2

line = repmat(point1,noPoints,1) + (repmat(vectorBase,noPoints,1) .* repmat((0:noPoints-1)',1,3));

end

