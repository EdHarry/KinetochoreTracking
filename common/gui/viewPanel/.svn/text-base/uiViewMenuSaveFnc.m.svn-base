function uiViewMenuSaveFnc
% callback function
axesH = get(gcbf,'CurrentAxes');
imgH = findobj(axesH,'Type','image','Parent',axesH);
data = get(imgH,'CData');

% Saves at 8 or 16bit? Info is stored in current object's UserData
par=get(gcbo,'UserData');
switch par
case '8bit'
   dpt=8;
case '16bit'
   dpt=16;
otherwise
end
   
% Select output file name
if(~isempty(data))
	[fName,dirName]=uiputfile('*.tif','View Panel: Save ...');
   if( isa(fName,'char') & isa(dirName,'char'))
      
      % The filename passed to imwrite must have the extension for imwrite to recognize the image format
      if isempty(findstr(fName,'.tif'))
         fName=strcat(fName,'.tif');
      end
      % Converts [0..1] data to 8 or 16 bit and then saves it to a file
      
      % check if image is 0..1 first
      if max(data(:)) > 1
          % make 0..1
          data = data - min(data(:));
          data = data ./max(data(:));
      end
      
      imwritend(data,[dirName,fName],dpt);
   end;
end;


