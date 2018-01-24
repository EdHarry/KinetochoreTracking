function area=cropStackNoCompression(area,path)
% cropStack lets the user crop a region of interest from an image stack
%
% SYNOPSIS   area=cropStack(area,path)
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

% Check input parameters
if nargin==0
   error('Not enough input arguments');
elseif nargin==1
   manualSave=1;
elseif nargin==2
    manualSave=0;
    if isempty(path) % In case the user passed an empty path instead of dropping it
        manualSave=1;
    end
else
    error('Wrong number of arguments');
end

% Creating output directory if needed
if manualSave==0 & exist(path)~=7
    % Directory does not exist - create it
    if path(2)==':' % Compatibility with Windows OS
        % Drive letter specified
        mkdir(path(1:3),path(4:end));
    else
        mkdir(path);
    end
    fprintf(1,'Directory %s successfully created.\n',path);
end

% Load First image
[fName,dirName] = uigetfile(...
    {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
    '*.tif','TIF files (*.tif)'
    '*.tiff','TIFF files (*.tiff)'
    '*.jpg;','JPG files (*.jpg)'
    '*.jpeg;','JPEG files (*.jpeg)'
    '*.*','All Files (*.*)'},...
    'Select first image');
if(isa(fName,'char') & isa(dirName,'char'))
    % Recover all file names from the stack
    outFileList=getFileStackNames([dirName fName]);
    % Number of files 
    n=length(outFileList);
else
    area=[];
    return
end

% The user can decide the number of images to be cropped
prompt={'Specify the number of images to be cropped'};
dlg_title='User input requested';
num_lines=1;
def={num2str(n)};
answer=fix(str2num(char(inputdlg(prompt,dlg_title,num_lines,def))));

% Check the selected number
if isempty(answer)
    area=[];
    return
end

if answer<1 | answer>n
    fprintf(1,'Invalid number of images specified. Using the default value (%d).\n',answer);
else
    % Crop outFileList
    n=answer;
    outFileList=outFileList(1:n);
end

% Read first image
imgOne=double(imread([dirName fName]));

if isempty(area)
    h=figure;
    set(h,'NumberTitle','off');
    set(h,'Name','Please draw region of interest');
    % Normalize imgOne
    imgOne=(imgOne-min(imgOne(:)))/(max(imgOne(:))-min(imgOne(:)));
    % Crop - if the user closes the window without drawing, roipoly will return an error
    try
        [imgCropped,area]=imcrop(imgOne);
    catch
        uiwait(msgbox('No polygon selected. Quitting','Error','modal'));
        area=[];
        return
    end
    % Close figure
    close(h);
    % Check selected polygon
    if area(3)==0 | area(4)==0
        uiwait(msgbox('Please chose an area, not a single pixel.','Error','error','modal'));
        area=[];
        return
    end
    % Round area (imcrop can give also non-integer boundaries)
    area=round(area);
    % Vertices
    y0=area(2); y=area(2)+area(4);
    x0=area(1); x=area(1)+area(3);
else
    % Vertices
    y0=area(1); y=area(3);
    x0=area(2); x=area(4);
end

% Check boundaries
if y0<=0, y0=1; end
if x0<=0, x0=1; end
if y>=size(imgOne,1), y=size(imgOne,1); end
if x>=size(imgOne,2), x=size(imgOne,2); end

% Select output directory
if manualSave==1
    path=uigetdir('','Select output directory');
    if path==0
        area=[];
        return
    end
end

% Initializing waitbar
h=waitbar(0,'Processing...');

% Processing files
for i=1:n
    
    % Current filename
    currentFile=char(outFileList(i));
    [fpath,fName,fno,fext]=getFilenameBody(currentFile);
   
    % Read image from disk
    img=imread(currentFile);
   
    % Cut image
    imgC=img(y0:y,x0:x);
   
    % Prepare filename with path
    filename=[path,filesep,'crop_',fName,fno,fext];
   
    % Write file to disk
    if strcmp(fext,'.tif') || strcmp(fext,'.tiff') || strcmp(fext,'.jpg') || strcmp(fext,'.jpeg')
        imwrite(imgC,filename, 'Compression', 'none');
    elseif strcmp(fext,'.bmp')
        imwrite(imgC,filename);     
    else
       disp('Unknown image format!!');
       return
    end
        
    % Update waitbar
    waitbar(i/n,h);
    
end

% Correct area
area=[y0 x0 y x];

% Close waitbar
close(h);
