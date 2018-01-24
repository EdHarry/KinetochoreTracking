function [x, y] = lsIntersectMaskLine(lines, mask_img, dir);

index = 1;
if dir == 1
   for i=1:length(lines)
      dI = find(diff(mask_img(:,lines(i))));
      if length(dI) > 0
         y(index:index+length(dI)-1) = dI;
         x(index:index+length(dI)-1) = lines(i);
         index = index+length(dI);
      end
   end
else
   for i=1:length(lines)
      dI = find(diff(mask_img(lines(i),:)));
      if length(dI) > 0
         x(index:index+length(dI)-1) = dI;
         y(index:index+length(dI)-1) = lines(i);
         index = index+length(dI);
      end
   end
end