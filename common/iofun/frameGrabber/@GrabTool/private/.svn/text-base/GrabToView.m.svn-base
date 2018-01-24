function gbT = GrabToView(gbT)
% grab first ...
[s,gbT.data] = mexFgInterface('grab');
leftV=squeeze(gbT.data(:,:,1))';
rightV=squeeze(gbT.data(:,:,2))';
% ... then show
if(~isempty(gbT.viewPanelH) & ishandle(gbT.viewPanelH))
   gbT.viewPanelH = uiViewPanelShowImg(leftV,0,gbT.viewPanelH);
else
   gbT.viewPanelH = uiViewPanelShowImg(leftV,1);
end;
if(~isempty(gbT.viewPanel2H) & ishandle(gbT.viewPanel2H))
   gbT.viewPanel2H = uiViewPanelShowImg(rightV,0,gbT.viewPanel2H);
else
   gbT.viewPanel2H= uiViewPanelShowImg(rightV,1);
end;

gbT.movie.data = [];
gbT.movie.cMap = [];
