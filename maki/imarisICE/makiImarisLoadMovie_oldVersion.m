function makiImarisLoadMovie_oldVersion( imarisHandle , dataStruct, decon )
%UNTITLED Loads a movie from a maki dataStruct (4d matrix saved in the
%.rawMoveName field) into imaris
% EHarry Nov 2011

if nargin < 3 || isempty(decon)
    decon = 1;
end

movie = dataStruct.rawMovieName;

if ~isnumeric(movie)
    movie = readOMEMatFile(fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName),[],[],decon);
end

movie = movie - min(movie(:));
movie = movie ./ max(movie(:));

if ~isfloat(movie) || size(size(movie),2) ~= 4
    error('Movie needs to be a 4d float to be loaded into Imaris');
end

aDataSet = imarisHandle.GetFactory.CreateDataSet;
aType = Imaris.tType.eTypeFloat;
aDataSet.Create(aType,size(movie,1),size(movie,2),...
    size(movie,3),1,size(movie,4));

% aDataSet.Create(aType,size(movie,1).*dataStruct.dataProperties.PIXELSIZE_XY,size(movie,2).*dataStruct.dataProperties.PIXELSIZE_XY,...
%     size(movie,3).*dataStruct.dataProperties.PIXELSIZE_Z,1,size(movie,4));


imarisHandle.SetDataSet(aDataSet);
imarisHandle.GetDataSet.SetExtendMinX(0);
imarisHandle.GetDataSet.SetExtendMinY(0);
imarisHandle.GetDataSet.SetExtendMinZ(0);
imarisHandle.GetDataSet.SetExtendMaxX(imarisHandle.GetDataSet.GetSizeX.*dataStruct.dataProperties.PIXELSIZE_XY);
imarisHandle.GetDataSet.SetExtendMaxY(imarisHandle.GetDataSet.GetSizeY.*dataStruct.dataProperties.PIXELSIZE_XY);
imarisHandle.GetDataSet.SetExtendMaxZ(imarisHandle.GetDataSet.GetSizeZ.*dataStruct.dataProperties.PIXELSIZE_Z);

imarisHandle.GetDataSet.SetChannelRange(0,0,1);

for i = 1:size(movie,4)
    temp = movie(:,:,:,i);
    temp = temp(:);
    imarisHandle.GetDataSet.SetDataVolumeAs1DArrayFloats(single(temp),0,i-1);
end


end

