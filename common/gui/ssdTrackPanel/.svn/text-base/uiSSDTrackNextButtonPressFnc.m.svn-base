function success = uiSSDTrackNextButtonPressFnc(pI,figH)
% callback for the NEXT button in the SSD Track panel
% the parameters pI and figH are optional; if they are given the function is 
% not working as a direct callback but as a part of a macro 
% started from another button (GO ON!) of the SSD Track panel

if(nargin == 2)
   ownH = findobj(gcbf,'Tag','UISSDTRACKNEXTBUTTON');
else
   figH = gcbf;
   ownH = gcbo;
   pI = get(gcbf,'UserData');
end;
goH = findobj(gcbf,'Tag','UISSDTRACKGOBUTTON');
editH = findobj(gcbf,'Tag','UISSDTRACKLOGEDIT');

% get the next position from the propagator
[pI.pos0, pI.shape0] = getPropagatorNextPos;
pI.posOff = [1,1];

if(isempty(pI.stack))
   [pI.sImg,pI.pos0,pI.posOff] = ...
      uiSSDTrackPanelInitSimg([],...
      pI.pos0,pI.posOff,pI.patchDim,pI.trackCmd);
   pI.nextIndx  = 1; % to keep the NEXT button enabled (see below)
   pI.nImgGrabbed = pI.nImgGrabbed + 1;
else
   [pI.sImg,pI.pos0,pI.posOff] = ...
      uiSSDTrackPanelInitSimg(pI.stack(:,:,pI.nextIndx),...
      pI.pos0,pI.posOff,pI.patchDim,pI.trackCmd);
end;

% start the matching
if(~isempty(pI.sImg))
   [pI.pos,dtls,dummy,newTimg] = ...
      ssdMatch(pI.tImg,pI.sImg,pI.tPos,pI.pos0,pI.shape0,pI.patchDim,pI.trackCmd);
   if(dtls.status<1)
      % exception handling for failed matching 
      pI = showResults(pI,pI.tPosRef,dtls,[],[]);
      msg=sprintf('Matching failed (%d),\n Manual initialization?',...
         dtls.status);
      ans = questdlg(msg,'SSD Tracking Failure','Yes','No','Yes');
      while( strcmp(ans,'Yes') & (dtls.status < 1))
         if(isempty(pI.stack))
            % since only minimal image information is grabbed online
            % we have to completely reinitialize a template
            [pI.tImg,pI.tPos,pI.tPosOff,pI.tPosRef,pI.patchDim,pI.trackCmd,pI.connectedView] = ...
               uiSSDTrackPanelInitTimg([],[],...
               pI.tImgMaskType,...
               pI.tImg,...
               pI.tPos,...
               pI.tPosRef,...
               pI.patchDim,...
               pI.trackCmd,...
               pI.connectedView,1);
            dtls.status = 0;
            % set the displacement orientation in the template coordinate frame 
            pI.tDispOri = pI.trackCmd.lori;
         else
            newPos0 = uiViewPanelSelectPt(pI.connectedView);
            [pI.sImg,pI.pos0,pI.posOff] = ...
               uiSSDTrackPanelInitSimg(pI.stack(:,:,pI.nextIndx),...
               newPos0,[1,1],pI.patchDim,pI.trackCmd);
            [pI.pos,dtls,dummy,newTimg] = ...
               ssdMatch(pI.tImg,pI.sImg,pI.tPos,pI.pos0,pI.shape0,...
               pI.patchDim,pI.trackCmd);
         end;
         if(dtls.status<0)
            pI = showResults(pI,pI.tPosRef,dtls,[],[]);
            msg=sprintf('Matching failed (%d),\n Manual initialization?',...
               dtls.status);
            ans = questdlg(msg,'SSD Track Failure','Yes','No','Yes');
         end;
      end;
   end;
   % update the structure for later matches
   if(dtls.status > 0)
      success = 1;
      oldtPosRef = pI.tPosRef;
      updateTimg = get(findobj(gcbf,'Tag','UISSDTRACKUPDATETIMGCHECK'),'Value');
      if(updateTimg == 1)
         % the new template image will be the filtered old search image
         pI.tImg = newTimg;
         pI.tPosOff = pI.posOff;
         pI.tPos = pI.pos;
         pI.trackCmd.lori = [];
         pI.tDispOri = [];
         pI.tPosRef = (dtls.shape*pI.tPosRef')';
         % update the template mask if there exists one
         if(pI.tImgMaskType > 1)
            % at the moment, this code is wrong  (GD Mai-3-1999):
            % 1.) It may happen that pI.nextIndx points to an inexistent frame
            %     (wrapped data loading is done later)
            % 2.) The mask should not be determined by a running through the 
            %     initialization again, but from propagation
            % => fix this code when switching to shape updated tracking
            error('"update template while a mask is set" is an invalid option in the current system; 5-3-99');
            %if(isempty(pI.stack))
            %   error('the option "update template while a mask is set" is not yet implemented; 10-22-98');
            %else
            %  [aux1,aux2,maskParams,pI.trackCmd.mask]=...
            %      uiSSDTrackPanelSetTmask(pI.stack(:,:,pI.nextIndx),...
            %      pI.posOff + pI.pos-1,pI.patchDim,pI.tImgMaskType);
            %   if(~isempty(maskParams))
            %      pI.tDispOri = maskParams(1);
            %      pI.trackCmd.lori = maskParams(1);
            %   end;
            %end;    
         end;   
      else
         if(isempty(pI.tDispOri))
            % there is not orientation stored with the template but may be some 
            % is returned by SSD match. If not dtls.ori = []
            % and pI.tDispOri remains []
            pI.tDispOri = dtls.lori;
         end;
         % the local orientation for the next track must be updated 
         pI.trackCmd.lori = rotLori(dtls.shape,pI.tDispOri);
      end; 
      % pI.pos0 = pI.pos;
      setPropagatorPosEstimate(pI.posOff+pI.pos-1,dtls.shape);
      % write the coordinates to the log file
      logFileName = get(editH,'String');
      if(~isempty(logFileName))
         if(isempty(pI.stack))
            uiSSDTrackPanelWriteLog([],dtls.status,...
               pI.posOff + pI.pos-1,oldtPosRef,...
               dtls.shape,dtls.lori,dtls.nPix,logFileName);
         else
            uiSSDTrackPanelWriteLog(pI.stackIndx(pI.nextIndx),...
               dtls.status,...
               pI.posOff + pI.pos-1,oldtPosRef,...
               dtls.shape,dtls.lori,dtls.nPix,logFileName);
         end;
      end;
      % show the position and template on screen      
      if(isempty(pI.stack))
         pI = showResults(pI,oldtPosRef,dtls,pI.trackCmd.lori,[]);
      else
         pI = showResults(pI,oldtPosRef,dtls,pI.trackCmd.lori,...
            pI.stackIndx(pI.nextIndx));
      end;
   else
      success = 0;
   end;
else
   success = 0;
end;

% set the text of the view panel
textH = findobj(pI.connectedView,'Type','uicontrol','Tag','UIVIEWTEXT');
if(isempty(pI.stack))
   set(textH,'String',sprintf('%d/%d',pI.nImgGrabbed,pI.nImgGrabbed));
else
   set(textH,'String',sprintf('%d/%d',pI.stackIndx(pI.nextIndx),...
      pI.stackIndx(1)-1+pI.stackTot));
end;
set(textH,'Visible','on');

% check the validity of the stack and set next index accordingly
if(~isempty(pI.stack))
   stackSze = size(pI.stack,3);
   if(pI.nextIndx < stackSze)
      % there are still images left in the stack
      pI.nextIndx = pI.nextIndx + 1;
   else
      % is the current index the absolutely last?
      if(pI.stackIndx(pI.nextIndx) == ...
            pI.stackIndx(1)-1+pI.stackTot)
         wrapBox = findobj(gcbf,'Tag','UISSDTRACKWRAPCHECK');
         wrap = get(wrapBox,'Value');
         if(wrap == 1)
            % only one wrap should be automatically made %
            set(wrapBox,'Value',0);
            % eventually load first part of the stack again
            if((pI.stackIndx(2)-1)~=pI.stackIndx(1))
               % check whether the next index is the image after the first
               % (might have been changed through the wrapping)
               [fpath,fbody,fno,fext] = getfilenamebody(pI.stackName);
               scndFileName = strcat(fpath,filesep,fbody,num2str(fno+1),fext);
               [auxStack,auxStackIndx,auxStackTot] = ...
                  imreadstack(scndFileName,...
                  uiSSDTrackPanelGetDflt('maxImgInStack')-1);
               pI.stack = cat(3,pI.stack(:,:,1),auxStack);
               pI.stackIndx = [pI.stackIndx(1),auxStackIndx];
            end;
            pI.nextIndx = 1;
         else
            % stop the tracking
            pI.nextIndx = [];
         end;  
      else
         % read the next part of the stack
         [fpath,fbody,fno,fext] = getfilenamebody(pI.stackName);
         nextFileName = strcat(fpath,filesep,fbody,...
            num2str(pI.stackIndx(pI.nextIndx)+1),fext);
         [auxStack,auxStackIndx,auxStackTot] = ...
            imreadstack(nextFileName,...
            uiSSDTrackPanelGetDflt('maxImgInStack')-1);
         pI.stack = cat(3,pI.stack(:,:,1),auxStack);
         pI.stackIndx = [pI.stackIndx(1),auxStackIndx];
         pI.nextIndx = 2;
      end;
   end;
end;

set(gcbf,'UserData',pI);
   
% en/disable the "next" / "go on!" button
if(~isempty(pI.nextIndx))
   set(ownH,'Enable','on');
   set(goH,'Enable','on');
else
   set(ownH,'Enable','off');
   set(goH,'Enable','off');
end;

%-------------------------------------------------------------
function newpI = showResults(pI,refPos,dtls,newOri,frmNo)

global ssdPlotLogFileName__;
newpI = pI;

if(~isempty(dtls.shape))
   shape = dtls.shape;
else
   shape = [1,0;0,1];
end;

rotRefPos = (shape*refPos')';

if(isempty(pI.stack))
   if(strcmp(pI.sImg.perm,'C'))
      pI.sImg = togglePermStatus(pI.sImg)
   end;
   newpI.connectedView = ...
      uiViewPanelShowImg(pI.sImg.data,0,pI.connectedView);
   winCtr = pI.pos;
else
   newpI.connectedView = ...
      uiViewPanelShowImg(pI.stack(:,:,pI.nextIndx),...
      0,pI.connectedView);
   winCtr = pI.pos + pI.posOff - 1;
end;
hold on;
cornerCrds = plotAffineRect(winCtr,...
   pI.patchDim,shape,'y-');
if(pI.tImgMaskType == 1)
   plot(winCtr(1)+rotRefPos(1),winCtr(2)+rotRefPos(2),'g+');
end;
plot(winCtr(1),winCtr(2),'g.');
if(~isempty(newOri))
   % the new orientation is already known because an orientation selective mask
   % was updated with the template update; or the template was not updated but 
   % the orientation was transformed to a new direction.
   plot([winCtr(1) + sqrt(sum(pI.patchDim).^2)/2*cos(newOri),...
         winCtr(1) - sqrt(sum(pI.patchDim).^2)/2*cos(newOri)],...
      [winCtr(2)+ sqrt(sum(pI.patchDim).^2)/2*sin(newOri),...
         winCtr(2) - sqrt(sum(pI.patchDim).^2)/2*sin(newOri)],'g-');
else
   % the new orientation is not yet known but there has been an orientation constraint applied.
   % Therefore show the old constraint as a dotted line 
   if(~isempty(dtls.lori))
   plot([winCtr(1) + sqrt(sum(pI.patchDim).^2)/2*cos(dtls.lori),...
         winCtr(1) - sqrt(sum(pI.patchDim).^2)/2*cos(dtls.lori)],...
      [winCtr(2)+ sqrt(sum(pI.patchDim).^2)/2*sin(dtls.lori),...
         winCtr(2) - sqrt(sum(pI.patchDim).^2)/2*sin(dtls.lori)],'g:');
   end;
end;
hold off;

if((dtls.status>0) & ~isempty(frmNo))
   % write the plot file
   fid = fopen(ssdPlotLogFileName__,'a');
   if(fid <0)
      return;
   end;
   
   % index of the frame
   fprintf(fid,'%6d ',frmNo);
   
   % coordinates of window center and reference point.
   fprintf(fid,'%8.2f %8.2f ',winCtr);
   fprintf(fid,'%8.2f %8.2f ',winCtr+rotRefPos);
   
   % x1 vector of the box coordinates
   fprintf(fid,'%8.2f %8.2f %8.2f %8.2f ',cornerCrds(1,:));
   % x2 vector of the box coordinates
   fprintf(fid,'%8.2f %8.2f %8.2f %8.2f ',cornerCrds(2,:));
   
   % direction vector
   if(~isempty(newOri))
      fprintf(fid,'%8.2f %8.2f ',...
         winCtr(1) + sqrt(sum(pI.patchDim).^2)/2*cos(newOri),...
         winCtr(1) - sqrt(sum(pI.patchDim).^2)/2*cos(newOri));
      fprintf(fid,'%8.2f %8.2f ',...
         winCtr(2)+ sqrt(sum(pI.patchDim).^2)/2*sin(newOri),...
         winCtr(2) - sqrt(sum(pI.patchDim).^2)/2*sin(newOri));
   else
      % no orientation is defined -> for column consitsincy fill 
      % in a zero vector arounf the winCtr
      fprintf(fid,'%8.2f %8.2f ',winCtr(1),winCtr(1));
      fprintf(fid,'%8.2f %8.2f ',winCtr(2),winCtr(2));
   end;
   
   fprintf(fid,'\n');
   fclose(fid);
end;

%-----------------------------------------------------------------------------
function newOri = rotLori(affMtx,oldOri)

newOri = oldOri;
if(isempty(oldOri))
   return;
end;

% the unit vector parallel to the old line is 
eLL = [sin(oldOri),-cos(oldOri)]';

% the transformed unit vector is
newELL = affMtx * eLL;
newELL = newELL / sqrt(sum(newELL.^2));

% rotate the new vector back by +90deg to be perpendicular to the line
newEPP = [-newELL(2),newELL(1)];

newOri = atan2(newEPP(2),newEPP(1));