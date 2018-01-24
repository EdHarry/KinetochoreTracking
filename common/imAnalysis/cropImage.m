function cropImage 
% cropImage allows the user to crop a region of interest from an image 
%
% SYNOPSIS   cropImage 
%
% INPUT      area [y0 x0 y x]
%                           y0 : y coordinate of the top-left corner
%                           x0 : x coordinate of the top-left corner
%                           y  : y coordinate of the bottom-right corner
%                           x  : x coordinate of the bottom-right corner
%            Pass area=[] to manually draw a region of interest
%
%            path (optional) output path - if path is not passed, the
%            user will be prompted to select the path through a save
%            dialog. If a directory is specified which does not exist, 
%            it will be created (if possible).
%
% AM Sept26 2003

% get current directory
currDir = cd;

enough = 0;
while ~enough   
    answer = questdlg('Crop an image?','','Yes','No','Yes');
    if strcmp (answer,'Yes')
        
        % Load image
        [fName,dirName] = uigetfile(...
            {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
            '*.tif','TIF files (*.tif)'
            '*.tiff','TIFF files (*.tiff)'
            '*.jpg;','JPG files (*.jpg)'
            '*.jpeg;','JPEG files (*.jpeg)'
            '*.*','All Files (*.*)'},...
            'Select image');
        
        % Change directory
        cd(dirName);
        
        % Read image to retrieve crop coordinates
        img=double(imread([dirName fName]));
        h=figure;
        set(h,'NumberTitle','off');
        set(h,'Name','Please draw region of interest');
        % Normalize img
        img=(img-min(img(:)))/(max(img(:))-min(img(:)));
        % Crop - if the user closes the window without drawing, roipoly will return an error
        try
            [imgCropped,area]=imcrop(img);
        catch
            uiwait(msgbox('No polygon selected. Quitting','Error','modal'));
            return
        end
        % Close figure
        close(h);
        % Round area (imcrop can give also non-integer boundaries)
        area=round(area);
        % Vertices
        y0=area(2); y=area(2)+area(4);
        x0=area(1); x=area(1)+area(3);
        
        % Check boundaries
        if y0<=0, y0=1; end
        if x0<=0, x0=1; end
        if y>=size(img,1), y=size(img,1); end
        if x>=size(img,2), x=size(img,2); end
        
        % Select output directory
        [file,path]=uiputfile('path','Select output directory');
        if path==0
            disp('Aborting...');
            return
        end
        
        % Read again image from disk
        img=imread([dirName fName]);
        
        % Cut image
        imgC=img(y0:y,x0:x);
        
        % Prepare filename with path
        filename=[path,filesep,'cut_',fName];
        
        % Write file to disk
        imwrite(imgC,filename);
        
    else 
        cd(currDir);
        enough = 1;
    end   
end % ~enough