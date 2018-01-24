function movie = readDecon(fileAndPath,metaData)
% EHarry Feb 2012

data = bfopenEdit(fileAndPath); % read file


% movie matrix goes [x,y,z,t,c]

[y,x] = size(data{1}{1,1});
sizeZ = str2double(data{4}.getImage(0).getPixels.getSizeZ.toString);

% movie = zeros(metaData.sizeX,metaData.sizeY,metaData.sizeZ,metaData.sizeT,metaData.sizeC);
movie = zeros(x,y,sizeZ,metaData.sizeT,metaData.sizeC);

if strcmp(metaData.dimOrder,'XYZCT') % only this dimORder for now
    count=1;
    for t = 1:metaData.sizeT
        for c = 1:metaData.sizeC
            for z = 1:sizeZ
                movie(:,:,z,t,c) = data{1}{count,1}'; % data goes (y,x) hence the transpose
                movie(:,:,z,t,c) = movie(:,end:-1:1,z,t,c);
                count=count+1;
            end
        end
    end
end


end

