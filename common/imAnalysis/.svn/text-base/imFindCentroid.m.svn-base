function centroid = imFindCentroid(x, y, radius, img)
% IMFINDCENTROID finds the centroid of a speckle 
%     
%
%
% SYNOPSIS      centroid = imFindCentroid(max, radius, img)
%
% INPUT                :         
% 
% OUTPUT               : 
%                           
% DEPENDENCES       prAlpha uses {    netAssemblyMaps
%                                     velocityMaps
%                                     prScoresInSeg
%                                       }
%
%                   prAlpha is used by { 
%                                           }
%
% Matthias Machacek 01/15/04
METHOD = 1;

[y_img, x_img] = size(img);

if METHOD == 1
    %use square patch
    
    if y-radius >=1 & y+radius <= y_img & x-radius >=1 & x+radius <= x_img
        %extract patch from image
        img_crop = img(y-radius:y+radius,x-radius:x+radius);
        s = size(img_crop);
        cx = 0;
        cy = 0;
        for l = 1 : s(2);
            cx = cx + sum(img(:,l).*l);
        end
        for l = 1 : s(1);
            cy = cy + sum(img(l,:).*l);
        end
        sTot=sum(img(:));
        centroid = [cx cy]/sTot;
        centroid(1) = x + centroid(1);
        centroid(2) = y + centroid(2);       
    else
        centroid = [x y];
    end
else 
    
    %use round patch
    for ix = max(2) - r : max(2) + r
        iy = round(sqrt(r^2 - (ix-max(2))^2) - max(2));
        int_sum = int_sum + img(iy, ix);
        xsum = xsum + ix * img(iy, ix);
        ysum = ysum + iy * img(iy, ix);   
    end
    
    centroid = 1;
    
end

