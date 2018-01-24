function uiSSDTrackPanelSet2PTOSCIparams()

global meanPixSize__;

params = getPropagatorParams;

prompt = {'Azimuth [deg]','Excursion [um]','Pixel Size [um]'};
title = 'SSD Track Oscillation Propagator';
lineNo = 1;

if(isempty(meanPixSize__))
   strVal3 = '';
else
   strVal3 = sprintf('%f',meanPixSize__);
end;

if(isempty(params))
   strVal1 = '';
   strVal2 = '';
else
   azi = 360 - params(1) * 180/pi;
   strVal1 = sprintf('%f',azi);
   if(~isempty(meanPixSize__))
      exc = params(2)*meanPixSize__;
      strVal2 = sprintf('%f',exc);
   else
      strVal2 = '';
   end;
end;

dfltVals = {strVal1,strVal2,strVal3};
result = inputdlg(prompt,title,lineNo,dfltVals);

if(isempty(result))
   return;
else
   azi = str2num(result{1});
   exc = str2num(result{2});
   pixSize = str2num(result{3});
   if(~isempty(pixSize))
      meanPixSize__ = pixSize;
      if(~isempty(exc))
         exc = exc / meanPixSize__;
      else
         return;
      end;
   else
      return;
   end;
   if(~isempty(azi))
      azi = azi/180*pi;
      if(azi < 0)
         azi = 2*pi + azi;
      end;
      azi = mod(azi,2*pi);
      azi = 2*pi - azi;
   else
      return;
   end;
   params = [azi,exc,0];
   setPropagatorParams(params);
end;

   




