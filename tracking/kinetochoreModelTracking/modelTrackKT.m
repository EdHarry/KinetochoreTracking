function [ output_args ] = modelTrackKT( dataStruct )
%MODELTRACKKT
%   EHarry Dec 2012

%turn warnings off
warningState = warning('off','all');

%get movie parameters
imageName = fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName);

% use deconvolved image or not
if isfield(dataStruct.dataProperties,'decon') && ~isempty(dataStruct.dataProperties.decon)
    decon = dataStruct.dataProperties.decon;
else
    decon = 0;
end

% get cropping if any
if isfield(dataStruct.dataProperties,'crop')
    crop = dataStruct.dataProperties.crop;
else
    crop = [];
end

%get initial guess of PSF sigma
param.psfSigma = dataStruct.dataProperties.psfSigma;

%get camera bit depth
if ~isfield(dataStruct.dataProperties,'bitDepth') || isempty(dataStruct.dataProperties.bitDepth)
    bitDepth = 16;
else
    bitDepth = dataStruct.dataProperties.bitDepth;
end

nFrames = dataStruct.dataProperties.movieSize(4);
initCoord = dataStruct.initCoord;
[~,tracksIdx] = convStruct2MatNoMS(dataStruct.tracks);
sisterList = dataStruct.sisterList(1).trackPairs(:,1:2);

movie = readOMEMatFile(imageName,1:nFrames,1,decon,crop);
[imageSizeX,imageSizeY,imageSizeZ,~] = size(movie);

% divide by bit depth
movie = movie ./ (2^bitDepth - 1);

% estimate bleaching
init = zeros(nFrames,1);
for iFrame = 1:nFrames
    amp = initCoord(iFrame).amp(:,1);
    init(iFrame) = robustMean(amp);
    clear amp
end
bleachingFactor = robustExponentialFit2(init);
bleachingFactor = bleachingFactor(2);
clear init

% loop frames
%for iFrame = 1:nFrames-1
for iFrame = 40
    % get target frame
    param.image = movie(:,:,:,iFrame+1);
    
    % initialise parameters
    baseCoord = [];
    amp = [];
    
    % initialise x0
    x0 = [];
    
    % find sisters in this frame
    pairs = [];
    for iSis = 1:size(sisterList,1)
        if tracksIdx(sisterList(iSis,1),iFrame) ~= 0 && tracksIdx(sisterList(iSis,2),iFrame) ~= 0 && tracksIdx(sisterList(iSis,1),iFrame+1) ~= 0 && tracksIdx(sisterList(iSis,2),iFrame+1) ~= 0
            pairs = [pairs; iSis];
            idx1_i = tracksIdx(sisterList(iSis,1),iFrame);
            idx2_i = tracksIdx(sisterList(iSis,2),iFrame);
            idx1_j = tracksIdx(sisterList(iSis,1),iFrame+1);
            idx2_j = tracksIdx(sisterList(iSis,2),iFrame+1);
            coord1_i = initCoord(iFrame).allCoordPix(idx1_i,1:3)';
            coord2_i = initCoord(iFrame).allCoordPix(idx2_i,1:3)';
            coord1_j = initCoord(iFrame+1).allCoordPix(idx1_j,1:3)';
            coord2_j = initCoord(iFrame+1).allCoordPix(idx2_j,1:3)';
            amp1 = initCoord(iFrame).amp(idx1_i,1);
            amp2 = initCoord(iFrame).amp(idx2_i,1);
            c1 = [coord1_i coord2_i];
            c2 = [coord1_j coord2_j];
            shift = mean(c2,2) - mean(c1,2);
            scale = sqrt(sum(diff(c2,1,2).^2)) / sqrt(sum(diff(c1,1,2).^2));
            c1N = (scale*eye(3))*(c1 - repmat(mean(c1,2),1,2));
            c2N = c2 - repmat(mean(c2,2),1,2);
            u = diff(c1N,1,2);
            v = diff(c2N,1,2);
            u = u / norm(u);
            v = v / norm(v);
            k = cross(u,v);
            costheta = dot(u,v);
            r = [ 0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
            r = costheta*eye(3) + r + k*k'*(1-costheta)/sum(k.^2);
            [psi,theta,phi] = eulerAnglesFromRotMat(r);
            % add to param
            baseCoord = cat(3,baseCoord,c1);
            amp = [amp; amp1 amp2];
            % add to x0
            x0 = [x0; shift; scale; phi; theta; psi];
        end
    end
    
    % add bleachingFactor
    x0(end+1) = exp(bleachingFactor);
    
    % add to param
    param.baseCoord = baseCoord;
    param.amp = amp;
    
    % calc a mean bg
    param.bg = mean(initCoord(iFrame).bg(:,1));
    
    % add no. of pairs
    param.numPairs = size(baseCoord,3);
    
    % fit
    [ solution, residuals, jacobian ] = modelTrackKT_fit( x0, param );
end

% turn wanings on
warning(warningState);

end

