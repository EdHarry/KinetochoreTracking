function setIDPropPosEstimate(pos, shape)
% service function for the ID and IDNOSHAPE propagator

global prop__;

prop__.currentPos = pos;
prop__.currentShape = shape;
return;
