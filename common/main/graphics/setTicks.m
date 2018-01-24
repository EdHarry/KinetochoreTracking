function setTicks(x,y,fmt)
%SETTICKS set the ticks and labels of the curent axis
%
% SYNOPSIS setTicks(x,y,fmt)
%
% INPUT x : vector with the ticks for the horizontal axis
%       y : vector with the ticks for the vertical axis
%       fmt: C-style printing format
%
% the change of the current ticks can be suppressed by entering []


if(~isempty(x))
   for i = 1:length(x)
      tl{i} = sprintf(fmt,x(i));
   end;
   set(gca,'XTick',x);
   set(gca,'XTickLabel',char(tl));
end;

if(~isempty(y))
   for i = 1:length(y)
      tl{i} = sprintf(fmt,y(i));
   end;
   set(gca,'YTick',y);
   set(gca,'YTickLabel',char(tl));
end;
