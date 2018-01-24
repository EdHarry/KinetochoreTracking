function st = changeVersion(fName,newV)
% change file version number
startSt = fName(1:regexp(fName,'_\d+\.mat'));
endSt = fName(regexp(fName,'.mat'):end);

st = [startSt int2str(newV) endSt];

end