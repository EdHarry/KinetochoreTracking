function gbT=SaveMpeg(gbT)

dirName = gbT.acqDir;
fileName = '*.mpg';
dfltName = strcat(dirName,filesep,fileName);

[fileName,dirName]=uiputfile(dfltName,...
   'Acquisition Panel: Save MPEG ...');

if( isa(fileName,'char') & isa(dirName,'char'))
   options(1) = gbT.movie.repeat;
   options(2) = 2;
   options(3) = 1;
   options(4) = 1;
   options(5) = 10;
   options(6) = 8;
   options(7) = 10;
   options(8) = 10;
   mpgwrite(gbT.movie.data,gbT.movie.cMap,...
      [dirName,fileName],options);
   gbT.acqDir = getFilenameBody([dirName,fileName]); 
end;


