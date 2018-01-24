function header = stk3dheader(rawMoviePath,rawMovieName)
%STK3DHEADER is a quick hack to read 3D+t STK files into the maki software package
%
%INPUT  rawMoviePath, rawMovieName: movie path and name.
%
%OUTPUT header: Same as output of readr3dheader.
%
%REMARKS: THIS IS A QUICK HACK - I HAVE HARDWIRED WHATEVER INFORMATION WAS
%NOT EASILY ACCESSIBLE
%
%Khuloud Jaqaman, August 2008

%read in first stack of time lapse
[stackData,numZSlices] = metaTiffRead(fullfile(rawMoviePath,rawMovieName));

%extract header information

%voxel size
header.pixelX = 0.129;
header.pixelY = 0.129;
header.pixelZ = 0.5;

%keep empty for now
header.firstImageAddress = []; %not sure what this is or what it is used for
header.lensID = 99999;

%movie size in X, Y and Z
header.numCols = stackData(1).width; %X in IMARIS
header.numRows = stackData(1).height; %Y in IMARIS
header.numZSlices = numZSlices;

%number of time points
fileList = searchFiles('.STK$',[],rawMoviePath,1);
header.numTimepoints = size(fileList,1);

%emission wavelength information
header.numWvs = 1;
header.zwtOrder = 'ztw';
header.wvl = 0.535;

%exposure time and neutral density filter
header.expTime = 0.1;
header.ndFilter = [];

%sampling time information
timeStamp = (0:header.numTimepoints-1)*5; %hardwire 5-sec sampling
timeStamp = repmat(timeStamp,numZSlices,1) + ...
    repmat((0:numZSlices-1)'*header.expTime,1,header.numTimepoints);
header.timestamp = timeStamp;
header.Time = timeStamp(:);

