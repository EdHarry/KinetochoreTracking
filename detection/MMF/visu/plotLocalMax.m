function plotLocalMax( movieParam, localMaxima , t )
%PLOTLOCALMAX plot locmax over images
% movieParam and localMaxima from movie analysis
% t time point
%   EHarry March 2012

% load frame
frame = readOMEMatFile(movieParam.imageName,t,movieParam.imageCh,movieParam.imageDecon,movieParam.imageCrop);

% get cand
cand = localMaxima(t).cands;

% get cord
cord = catStruct(1,'cand.Lmax');

% draw max intensity projection of the frame
imshow(max(frame,[],3),[]);
hold on

% plot coords
plot(cord(:,2),cord(:,1),'g*');

end

