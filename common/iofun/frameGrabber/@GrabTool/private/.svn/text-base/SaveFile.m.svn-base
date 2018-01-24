function gbT = SaveFile(gbT)

newDirString = [];
newFileString = [];

dirName = gbT.acqDir;
fileName = gbT.acqFile;
if(isempty(fileName))
   fileName = '*.tif';
end;
dfltName = strcat(dirName,filesep,fileName);

% save
if(~isempty(gbT.data))
   if(~isempty(gbT.movie.stackIndx))
      [fileName,dirName]=uiputfile(dfltName,...
         'Acquisition Panel: Save Stack ...');
      
      % The filename passed to imwrite must have the extension for imwrite to recognize the image format
      if isempty(findstr(fileName,'.tif'))
         fileName=strcat(fileName,'.tif');
      end

      if( isa(fileName,'char') & isa(dirName,'char'))
			% Select bit depth
         choice=Questdlg('Select bit depth...','Input requested','8 bit','16 bit','8 bit');
         switch choice
         case '8 bit'
            dpt=8;
         case '16 bit'
            dpt=16;
         otherwise
            error('You should choose between 8 and 16 bit');
         end
         % Writes
         imwritestacknd(gbT.data,[dirName,fileName],dpt);
         [newDirString,fbody,fno,fext] = ...
            getFilenameBody([dirName,fileName]);
         if(~isempty(fno))
            newFileString = ...
               strcat(fbody,int2str(str2num(fno)+length(gbT.movie.stackIndx)), ...
               fext);
         end;
   	end;      
   else
      [fileName,dirName]=uiputfile(dfltName,...
         'Acquisition Panel: Save ...');
      
      % The filename passed to imwrite must have the extension for imwrite to recognize the image format
      if isempty(findstr(fileName,'.tif'))
         fileName=strcat(fileName,'.tif');
      end
      
      if( isa(fileName,'char') & isa(dirName,'char'))
         % Select bit depth
         choice=Questdlg('Select bit depth...','Input requested','8 bit','16 bit','8 bit');
         switch choice
         case '8 bit'
            dpt=8;
         case '16 bit'
            dpt=16;
         otherwise
            error('You should choose between 8 and 16 bit');
         end
         
         % Writes
         imwritend(gbT.data,[dirName,fileName],dpt);
         [newDirString,fbody,fno,fext] = ...
            getFilenameBody([dirName,fileName]);
         if(~isempty(fno))
            newFileString = strcat(fbody,int2str(str2num(fno)+1),fext);
         else
            newFileString = strcat(fbody,'1',fext);
         end;
      end;
   end;
end;

% update the directory and file edit fields
if( isa(newFileString,'char') & isa(newDirString,'char'))
   gbT.acqDir=newDirString;
   gbT.acqFile=newFileString;
end;
