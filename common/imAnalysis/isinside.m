function ans = isinside(smallRect,rect,clip)
%ISINSIDE checks whether one rectangular area is completely inside another
%
% SYNOPSIS ans = isinside(smallRect,rect,clip)
%
% INPUT  smallRect : defintion of rectangular area which should be inside of
%        rect :
%               both rectangular areas are defined as 
%               [ul(1), ul(2), width, height]
%        clip : (optional) default 0
%               1 -> smallRect is clipped such that it is completely
%                    inside; the result is written to out
%               0 -> just test
% 
% OUTPUT ans : dependent on clip either a new rectangular area
%              or a boolean answer
%
if(nargin < 3)
   clip = 0;
end;

smallRectLr(1) = smallRect(1) + smallRect(3);
smallRectLr(2) = smallRect(2) + smallRect(4);
rectLr(1) = rect(1) + rect(3);
rectLr(2) = rect(2) + rect(4);

ans = smallRect(1)>=rect(1);
ans = ans & (smallRect(2)>=rect(2));
ans = ans & (smallRectLr(1) <=rectLr(1));
ans = ans & (smallRectLr(2) <=rectLr(2));

if(clip)
   outside = ~ans;
   ans = smallRect;
   if(outside)
      if(ans(1) < rect(1)) ans(1) = rect(1); end;
      if(ans(2) < rect(2)) ans(2) = rect(2); end;
      if(smallRectLr(1) > rectLr(1)) 
         ans(3) = rectLr(1) - ans(1); 
      end;
      if(smallRectLr(2) > rectLr(2)) 
         ans(3) = rectLr(2) - ans(2);
      end;
   end;
end;




