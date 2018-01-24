function makiImarisLoadMovie( imarisHandle , dataStruct, decon,c2Load )
%UNTITLED Loads a movie from a maki dataStruct (4d matrix saved in the
%.rawMoveName field) into imaris
% EHarry Nov 2011

if nargin < 3 || isempty(decon)
    decon = 1;
end

if nargin < 4
    c2Load = [];
end

movie = dataStruct.rawMovieName;

if ~isnumeric(movie)
    movie = readOMEMatFile(fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName),[],c2Load,decon);
end

noC = size(movie,5);

chMinMax = zeros(noC,2);
for c = 1:noC
    temp = movie(:,:,:,:,c);
    chMinMax(c,1) = min(temp(:));
    chMinMax(c,2) = max(temp(:));
    %temp = temp - min(temp(:));
    %temp = temp ./ max(temp(:));
    %movie(:,:,:,:,c) = temp;
end

if ~isfloat(movie) || (size(size(movie),2) ~= 4 && size(size(movie),2) ~= 5 && size(size(movie),2) ~= 3)
    error('Movie needs to be a 3d, 4d or 5d float to be loaded into Imaris');
end

aDataSet = imarisHandle.GetFactory.CreateDataSet;
aType = Imaris.tType.eTypeFloat;
aDataSet.Create(aType,size(movie,1),size(movie,2),...
    size(movie,3),size(movie,5),size(movie,4));

% aDataSet.Create(aType,size(movie,1).*dataStruct.dataProperties.PIXELSIZE_XY,size(movie,2).*dataStruct.dataProperties.PIXELSIZE_XY,...
%     size(movie,3).*dataStruct.dataProperties.PIXELSIZE_Z,1,size(movie,4));


imarisHandle.SetDataSet(aDataSet);
imarisHandle.GetDataSet.SetExtendMinX(0);
imarisHandle.GetDataSet.SetExtendMinY(0);
imarisHandle.GetDataSet.SetExtendMinZ(0);
imarisHandle.GetDataSet.SetExtendMaxX(imarisHandle.GetDataSet.GetSizeX.*dataStruct.dataProperties.PIXELSIZE_XY);
imarisHandle.GetDataSet.SetExtendMaxY(imarisHandle.GetDataSet.GetSizeY.*dataStruct.dataProperties.PIXELSIZE_XY);
imarisHandle.GetDataSet.SetExtendMaxZ(imarisHandle.GetDataSet.GetSizeZ.*dataStruct.dataProperties.PIXELSIZE_Z);


for c = 1:noC
    imarisHandle.GetDataSet.SetChannelRange(c-1,chMinMax(c,1),chMinMax(c,2));
end

for j = 1:size(movie,5)
    for i = 1:size(movie,4)
        temp = movie(:,:,:,i,j);
        temp = temp(:);
        imarisHandle.GetDataSet.SetDataVolumeAs1DArrayFloats(single(temp),j-1,i-1);
    end
end

end

