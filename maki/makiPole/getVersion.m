function v = getVersion(fName)
% Service function to retrieve the current version from a filename in the
% job data structure
%
% The function assumes that the version index is the last number before
% '.mat' and is preceeded by an '_';
%
% for example: 'gaga_cpi11_anythingelse_23.mat'
% will generate a numerical value 23
v = str2num(fName(regexp(fName,'_\d+\.mat')+1:regexp(fName,'.mat')-1));

end