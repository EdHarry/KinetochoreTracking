function plotImageWithFeatures(image,positions)

%make 3 layers out of normalized image
image = double(image);
imageNorm = image/max(image(:));
imageN3 = repmat(imageNorm,[1 1 3]);

%get number of pixel in each side of image
[numPixelsX,numPixelsY] = size(image);

%get set of colors for labeling
color = [1 0 0];

%label maxima
for i=1:size(positions,1)
    pos = (round(positions(:,1)-1))*numPixelsX + round(positions(:,2));
    %     pos = [pos pos+1 pos+numPixelsX pos+1+numPixelsX];
    for j=1:3
        imageN3(pos+(j-1)*numPixelsX*numPixelsY)=color(j);
    end
end

%plot image
imtool(imageN3);
