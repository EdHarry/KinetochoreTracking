function [tImg,tPos,tPosOff,tPosRef,tDim,newTrackCmd,newViewH] = ...
   uiSSDTrackPanelInitTimg(img,idx,maskType,oldTimg,oldTpos,oldTposRef,oldTdim,trackCmd,viewH,silent)
% service function for the interactive initialization of the 
% template
%
% the input parameter silent is optional; if 1 init template function 
% automatically uses the old template (if available) to initialize a
% new template image
global ssdPlotLogFileName__;
global bwOverlays__;
global debuggingMode__;

automatic = 0;

% tPosRef is in principle the center of the template unless we compute 
% center with a mask

tPosRef = [0,0];

if(nargin<10)
   silent = 0;
end;

% restore a copy of the trackCmd in newTrackCmd
newTrackCmd = trackCmd;
newTrackCmd.lori = [];   % remove possibly existing old local orientation 
newTrackCmd.mask = [];   % remove possibly existing template mask 

% get the image
if(isempty(img) & isReadyFg)
   img = grab;
   idx = [];
else
   if(isempty(img) & ~isReadyFg)
      return;
   end;
end;

% show the image
newViewH=uiViewPanelShowImg(img,1,viewH);
if(newViewH < 0) return; end;

if(isempty(oldTimg))
   % there is no old template which could be used for initialization purposes
   rectSelect = 1;
else
   % there is an old template which might be used to reinitialize
   if(silent == 1)
      rectSelect = 0; 
   else	
      ans = questdlg('Use existing template?',...
         'SSD Track Panel Init','Yes','No','Yes');
      if(strcmp(ans,'Yes'))
         rectSelect = 0;
      else
         rectSelect = 1;
      end;
   end;
end;

if(rectSelect == 1)
   % select rect interactively
   sRect = [0,0,0,0];
   while( (sRect(3) < 2) | (sRect(4) < 2))
      sRect = uiViewPanelSelectRect(newViewH);
   end;
   sRect = round(sRect);
   tDim = 2*floor(sRect(3:4)/2);
   tPos = sRect(1:2) + tDim/2;
   if(maskType > 0)
      [tPos,tPosRef,maskParams,newTrackCmd.mask]=...
         uiSSDTrackPanelSetTmask(img,tPos,tDim,maskType,newViewH);
      if(~isempty(maskParams))
         newTrackCmd.lori = maskParams(1);  
      end;
   end;
else
   % try to do SSD matching between old template and new template image
   auxPos0 = uiViewPanelSelectPt(newViewH);
   [auxImg,auxPos0,auxPosOff] = ...
      uiSSDTrackPanelInitSimg(img,...
      auxPos0,[1,1],oldTdim,trackCmd);
   if(~isempty(auxImg))
      trackCmd.lori = []; % remove possibly existing old local orientation, 
                          % which is obsolete to initialize new template
      [auxPos,dtls] = ...
         ssdMatch(oldTimg,auxImg,oldTpos,auxPos0,[],oldTdim,trackCmd);
      if(dtls.status<1)
         % attempt failed ...
         msg=sprintf('Matching failed (%d),\n switching to manual initialization',...
            dtls.status);
         msgbox(msg,'SSD Track Message','modal');
         sRect = [0,0,0,0];
         while( (sRect(3) < 2) | (sRect(4) < 2))
            sRect = uiViewPanelSelectRect(newViewH);
         end;
         sRect = round(sRect);
         tDim = 2*floor(sRect(3:4)/2);
         tPos = sRect(1:2) + tDim/2;
         if(maskType > 0)
            [tPos,tPosRef,maskParams,newTrackCmd.mask]=...
               uiSSDTrackPanelSetTmask(img,tPos,tDim,maskType,newViewH);
            if(~isempty(maskParams))
               newTrackCmd.lori = maskParams(1);
            end;  
         end;
      else
         % matching was successful
         tPos = auxPosOff + auxPos - 1;
         tPosRef = (dtls.shape*oldTposRef')';
         tDim = oldTdim;
         automatic = 1;
         if(~isempty(newTrackCmd.mask))
            % update mask at new template position
            [auxPos,auxPosRef,maskParams,newTrackCmd.mask]=...
               uiSSDTrackPanelSetTmask(img,tPos,tDim,maskType,newViewH);
            if(~isempty(maskParams))
               newTrackCmd.lori = maskParams(1);
            end;            
         end;
      end;
   end;
end;

figure(newViewH);
hold on;

if(~isempty(newTrackCmd.mask))
   if(~bwOverlays__)
      plotmask(newTrackCmd.mask,tPos,'m.');
   else
      plotmask(newTrackCmd.mask,tPos,'k.');
   end;
end;

if(~isempty(newTrackCmd.lori))
   if(~bwOverlays__)
      plot([tPos(1) + sqrt(sum(tDim).^2)/2*cos(newTrackCmd.lori),...
            tPos(1) - sqrt(sum(tDim).^2)/2*cos(newTrackCmd.lori)],...
         [tPos(2)+ sqrt(sum(tDim).^2)/2*sin(newTrackCmd.lori),...
            tPos(2)-sqrt(sum(tDim).^2)/2*sin(newTrackCmd.lori)],'g-');
   else
      plot([tPos(1) + sqrt(sum(tDim).^2)/2*cos(newTrackCmd.lori),...
            tPos(1) - sqrt(sum(tDim).^2)/2*cos(newTrackCmd.lori)],...
         [tPos(2)+ sqrt(sum(tDim).^2)/2*sin(newTrackCmd.lori),...
            tPos(2)-sqrt(sum(tDim).^2)/2*sin(newTrackCmd.lori)],'k-');   
   end;
end;

switch maskType,
case 0, plot(tPos(1)+tPosRef(1),tPos(2)+tPosRef(2),'b+');
   plot(tPos(1),tPos(2),'b.');     
case 1, plot(tPos(1)+tPosRef(1),tPos(2)+tPosRef(2),'b+');
   plot(tPos(1),tPos(2),'b.');
otherwise, 
   if(~bwOverlays__)
%      plot(tPos(1)+tPosRef(1),tPos(2)+tPosRef(2),'r+');
      plot(tPos(1),tPos(2),'g.');
   else
%     plot(tPos(1)+tPosRef(1),tPos(2)+tPosRef(2),'k+');
      plot(tPos(1),tPos(2),'w.');
   end;
end;

cDim = tDim + ...
   2*ceil(uiSSDTrackPanelGetDflt('bdFactor')*trackCmd.sigma);
cRect = [round(tPos)-cDim/2,cDim];

if(~bwOverlays__)
   cornerCrds = plotrect([tPos-tDim/2,tDim],'y-');
   if(debuggingMode__)
      plotrect(cRect,'y-');
   end;
else
   cornerCrds = plotrect([tPos-tDim/2,tDim],'k-');
end;

hold off;

if(~isempty(idx))
   % write the plot coordinates out to the plot file
   fid = fopen(ssdPlotLogFileName__,'a');   
   
   % image index
   fprintf(fid,'%6d ',idx);

   % coordinates of box center and reference point.
   fprintf(fid,'%8.2f %8.2f ',tPos);
   fprintf(fid,'%8.2f %8.2f ',tPos+tPosRef);
   
   % x1 vector of the box coordinates
   fprintf(fid,'%8.2f %8.2f %8.2f %8.2f ',cornerCrds(1,:));
   % x2 vector of the box coordinates
   fprintf(fid,'%8.2f %8.2f %8.2f %8.2f ',cornerCrds(2,:));
   
   % direction vector
   if(~isempty(newTrackCmd.lori))
      fprintf(fid,'%8.2f %8.2f ',...
         tPos(1) + sqrt(sum(tDim).^2)/2*cos(newTrackCmd.lori),...
         tPos(1) - sqrt(sum(tDim).^2)/2*cos(newTrackCmd.lori));
      fprintf(fid,'%8.2f %8.2f ',...
         tPos(2)+ sqrt(sum(tDim).^2)/2*sin(newTrackCmd.lori),...
         tPos(2) - sqrt(sum(tDim).^2)/2*sin(newTrackCmd.lori));
   else
      % no orientation is defined -> for column consitsincy fill 
      % in a zero vector arounf the winCtr
      fprintf(fid,'%8.2f %8.2f ',tPos(1),tPos(1));
      fprintf(fid,'%8.2f %8.2f ',tPos(2),tPos(2));
   end;
   
   fprintf(fid,'\n');
   fclose(fid);
end;
   
% crop the image
tImg.data = imcrop(img,cRect);
tImgSze = size(tImg.data);
if( ((tImgSze(1) - 1)~=cRect(4)) | ((tImgSze(2) - 1)~=cRect(3)))
   tImg.data = [];
else
   tImg.perm = 'M';
   tImg.prefilter = 1;
end;

% set the coordinates
tPos = (tPos - cRect(1:2)) + 1;
tPosOff = cRect(1:2);

nextButtonH = findobj(gcbf,'Tag','UISSDTRACKNEXTBUTTON');
goButtonH = findobj(gcbf,'Tag','UISSDTRACKGOBUTTON');
editH = findobj(gcbf,'Tag','UISSDTRACKLOGEDIT');
if(~isempty(tImg.data))
   % enable the "next" / "go on!" buttons
   set(nextButtonH,'Enable','on');
   set(goButtonH,'Enable','on');  
   % write the coordinates to the log file
   logFileName = get(editH,'String');
   if(~isempty(logFileName))
      if(automatic == 0)
         uiSSDTrackPanelWriteLog(idx,0,...
            tPosOff + tPos-1,tPosRef,...
            [1,0;0,1],newTrackCmd.lori,0,logFileName);
      else
         uiSSDTrackPanelWriteLog(idx,-dtls.status,...
            tPosOff + tPos-1,tPosRef,...
            dtls.shape,newTrackCmd.lori,dtls.nPix,logFileName);
      end;
   end;
else
   set(nextButtonH,'Enable','off');
end;


