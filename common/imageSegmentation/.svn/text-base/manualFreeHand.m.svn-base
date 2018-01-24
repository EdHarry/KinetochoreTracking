function varargout = manualFreeHand(varargin)
% MANUALFREEHAND M-file for manualFreeHand.fig
%      MANUALFREEHAND is a GUI for manual segmentation
%
%      MANUALFREEHAND, by itself, creates a new MANUALFREEHAND or raises the existing
%      singleton*.
%
%      H = MANUALFREEHAND returns the handle to a new MANUALFREEHAND or the handle to
%      the existing singleton*.
%
%      MANUALFREEHAND('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANUALFREEHAND.M with the given input arguments.
%
%      MANUALFREEHAND('Property','Value',...) creates a new MANUALFREEHAND or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before manualSegPanel_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to manualFreeHand_OpeningFcn via
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help manualFreeHand

% Last Modified by GUIDE v2.5 08-Oct-2008 14:12:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @manualFreeHand_OpeningFcn, ...
                   'gui_OutputFcn',  @manualFreeHand_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --------------------------------------------------------------------
function menu_close_Callback(hObject, eventdata, handles)
% hObject    handle to menu_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
user_response = yesnodialog('Title','Confirm Close!');
switch user_response
case {'No'}
    % take no action
case 'Yes'
    % Prepare to close GUI application window
    %                  .
    %                  .
    %                  .
    delete(handles.figure1)
end

% --- Executes just before manualFreeHand is made visible.
function manualFreeHand_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to manualFreeHand (see VARARGIN)

handles.directory_name = pwd;

maskinfo = update_mask_num(hObject, handles);
handles.maskinfo = maskinfo;

val = get(handles.menu_color,'Value');
switch val
case 1
    handles.colormap_type = 'gray';
case 2
    handles.colormap_type = 'Jet';
case 3
    handles.colormap_type = 'HSV';
case 4
    handles.colormap_type = 'HOT';    
end

% Choose default command line output for manualFreeHand
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes manualFreeHand wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = manualFreeHand_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ------------------------------------------------------------
% Update number of masks under the handles.directory_name
% ------------------------------------------------------------
function maskinfo = update_mask_num(hObject, handles)

FileType = '*.tif';

fullmaskpath = [handles.directory_name filesep FileType];
maskdir = dir(fullmaskpath);
mask_num = length(maskdir);

maskinfo.fullmaskpath = fullmaskpath;
maskinfo.maskdir = maskdir;
maskinfo.mask_num = mask_num;

set(handles.imgnumtext, 'String', mask_num);

if(mask_num > 0)
    set(handles.loadbutton1, 'Enable', 'on'); 
    set(handles.slider1, 'Enable', 'on');    
    set(handles.savebutton, 'Enable', 'on');    
    set(handles.startbutton, 'Enable', 'on');    
else
    set(handles.loadbutton1, 'Enable', 'off');
    set(handles.slider1, 'Enable', 'off');    
    set(handles.savebutton, 'Enable', 'off');
    set(handles.startbutton, 'Enable', 'off');    
end


% handles.fullmaskpath = fullmaskpath;
% handles.maskdir = maskdir;
% handles.mask_num = mask_num;
% set(handles.imgnumtext, 'String', mask_num);
% 
% if(handles.mask_num > 0)
%     set(handles.loadbutton1, 'Enable', 'on');
% else
%     set(handles.loadbutton1, 'Enable', 'off');    
% end
% 
% % Update handles structure
% guidata(hObject, handles);


% --- Executes on button press in browsebutton1.
function browsebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to browsebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


directory_name = uigetdir;
if directory_name == 0
    return;
end  
handles.directory_name = directory_name;
handles.result_directory_name = directory_name;

set(handles.inputdirectoryedit,'String',handles.directory_name);
maskinfo = update_mask_num(hObject, handles);
handles.maskinfo = maskinfo;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inputdirectoryedit_Callback(hObject, eventdata, handles)
% hObject    handle to inputdirectoryedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputdirectoryedit as text
%        str2double(get(hObject,'String')) returns contents of inputdirectoryedit as a double
new_directory_name = get(hObject,'String');

if isdir(new_directory_name)
    handles.directory_name = new_directory_name;
    handles.result_directory_name = new_directory_name;
    set(hObject,'String', handles.directory_name);      
    maskinfo = update_mask_num(hObject, handles);
    handles.maskinfo = maskinfo;
    guidata(hObject, handles);
else
    errordlg([new_directory_name ' is not a valid directory name'],'Bad Input','modal')
    set(hObject,'String', handles.directory_name);  
end


% --- Executes during object creation, after setting all properties.
function inputdirectoryedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputdirectoryedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',pwd);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
mask_idx = get(handles.slider1,'Value');
mask_idx = round(mask_idx);
set(handles.slider1,'Value',mask_idx);
thismask = imread([handles.directory_name filesep handles.maskinfo.maskdir(mask_idx).name]);
axes(handles.axes1);
imshow(thismask, []);   colormap(handles.colormap_type);
text(handles.colSize+5,20,['Image Index:' num2str(mask_idx)],'BackgroundColor',[.7 .9 .7]);

handles.thismask = thismask;
handles.mask_idx = mask_idx;
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in loadbutton1.
function loadbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FileType = '*.tif';
MNfullmaskpath = [handles.result_directory_name filesep 'manual_masks' filesep FileType];

MNmaskdir = dir(MNfullmaskpath);
MNmask_num = length(MNmaskdir);
i = 1;
j = 1;
flag = 0;
mask_idx = 1;
while i<=MNmask_num & j <= handles.maskinfo.mask_num
    if strcmp( MNmaskdir(i).name,handles.maskinfo.maskdir(j).name ) == 1
        i = i + 1;
        j = j + 1;    
    else
        mask_idx = i;
        break;
    end
    if j == handles.maskinfo.mask_num;
    	flag = 1;
        mask_idx = 1;        
    end
end
slider_step(1) = 1/handles.maskinfo.mask_num;
slider_step(2) = 1/handles.maskinfo.mask_num;
thismask = imread([handles.directory_name filesep handles.maskinfo.maskdir(mask_idx).name]);
set(handles.slider1, 'Visible', 'on','Max',handles.maskinfo.mask_num,...
    'Value',1,'Min',1,'sliderstep',slider_step);

axes(handles.axes1);    imshow(thismask, []);   colormap(handles.colormap_type);

set(handles.startbutton, 'Visible', 'on', 'Enable', 'on');
set(handles.savebutton, 'Visible', 'on', 'Enable', 'on');

handles.thismask = thismask;
handles.mask_idx = mask_idx;

handles.rowSize = size(thismask,1);
handles.colSize = size(thismask,2);

text(handles.colSize+5,20,['Image Index:' num2str(mask_idx)],'BackgroundColor',[.7 .9 .7]);

% Update handles structure
guidata(hObject, handles);

if flag == 1;
    set(handles.statusbar, 'String', ['Status: All the image has been segmented.  Please run protrusion analysis']);    
    set(handles.slider1, 'Enable', 'off');    
    set(handles.savebutton, 'Enable', 'off');
    set(handles.startbutton, 'Enable', 'off');        
else
    set(handles.statusbar, 'String', ['Status: Images are loaded for manual segmentation']);
end


% --- Executes during object creation, after setting all properties.
function banneraxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to banneraxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate banneraxes
set(hObject, 'Units', 'pixels');
imagePath = ['banner.jpg'];

%imagePath = ['Chou.jpg'];
handles.banner = imread(imagePath);
[info.Height info.Width info.Channel] = size(handles.banner);
position = get(hObject, 'Position');
axes(hObject);
image(handles.banner)
set(hObject, ...
    'Visible', 'off', ...
    'Units', 'pixels', ...
    'Position', [13 450 info.Width info.Height]);


% --------------------------------------------------------------------
function menu_Callback(hObject, eventdata, handles)
% hObject    handle to menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

output_dir = [handles.result_directory_name filesep 'manual_masks'];
savepath = [output_dir filesep handles.maskinfo.maskdir(handles.mask_idx).name];
imwrite(handles.segimg,savepath,'tif');
set(handles.statusbar, 'String', ['Status: Mask ' savepath ' has been saved']);
pause(1);
if handles.mask_idx < handles.maskinfo.mask_num
    handles.mask_idx = handles.mask_idx + 1;
    handles.thismask = imread([handles.directory_name filesep handles.maskinfo.maskdir(handles.mask_idx).name]);

    set(handles.slider1,'Value',handles.mask_idx);    
    
else
    set(handles.statusbar, 'String', ['Status: All the image has been segmented.  Please run protrusion analysis']);    
    set(handles.slider1, 'Enable', 'off');    
    set(handles.savebutton, 'Enable', 'off');
    set(handles.startbutton, 'Enable', 'off');            
    return;
end

guidata(hObject, handles);
startbutton_Callback(hObject, eventdata, handles);

% --- Executes on button press in startbutton.
function startbutton_Callback(hObject, eventdata, handles)
% hObject    handle to startbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.startbutton, 'Visible', 'off');
set(handles.savebutton, 'Visible', 'off');

set(handles.statusbar, 'String', ['Status: Segment Image ' num2str(handles.mask_idx)]);

output_dir = [handles.result_directory_name filesep 'manual_masks'];

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end


axes(handles.axes1);
img = handles.thismask;
imshow(img, []);    colormap(handles.colormap_type);
text(handles.colSize+5,20,['Image Index:' num2str(handles.mask_idx)],'BackgroundColor',[.7 .9 .7]);

% Update handles structure
guidata(hObject, handles);
img_show(:,:,1) = zeros(size(img));
img_show(:,:,2) = mat2gray(img);
img_show(:,:,3) = mat2gray(img);
h = imfreehand;
position = wait(h);     
set(handles.startbutton, 'Visible', 'on');
set(handles.savebutton, 'Visible', 'on');

valid_position = position;
valid_position(find(position(:,1)<0.5),1) = 1;
valid_position(find(position(:,2)<0.5),2) = 1;
valid_position(find(position(:,1)>handles.colSize),1) = handles.colSize;
valid_position(find(position(:,2)>handles.rowSize),2) = handles.rowSize;
segimg = roipoly(img,valid_position(:,1),valid_position(:,2));

handles.segimg = segimg;
img_show(:,:,1) = 0.4*mat2gray(segimg);
imshow(img_show);   colormap(handles.colormap_type);
text(handles.colSize+5,20,['Image Index:' num2str(handles.mask_idx)],'BackgroundColor',[.7 .9 .7]);

guidata(hObject, handles);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes on selection change in menu_color.
function menu_color_Callback(hObject, eventdata, handles)
% hObject    handle to menu_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menu_color contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_color

val = get(handles.menu_color,'Value');
switch val
case 1
    handles.colormap_type = 'gray';
case 2
    handles.colormap_type = 'Jet';
case 3
    handles.colormap_type = 'HSV';
case 4
    handles.colormap_type = 'HOT';    
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function menu_color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


