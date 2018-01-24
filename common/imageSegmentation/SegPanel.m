function varargout = SegPanel(varargin)
% SegPanel M-file for SegPanel.fig
%      SegPanel, by itself, creates a new SegPanel or raises the existing
%      singleton*.
%
%      H = SegPanel returns the handle to a new SegPanel or the handle to
%      the existing singleton*.
%
%      SegPanel('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SegPanel.M with the given input arguments.
%
%      SegPanel('Property','Value',...) creates a new SegPanel or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before The_Crapper_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SegPanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SegPanel

% Last Modified by GUIDE v2.5 06-May-2009 12:54:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SegPanel_OpeningFcn, ...
                   'gui_OutputFcn',  @SegPanel_OutputFcn, ...
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

% --- Executes just before SegPanel is made visible.
function SegPanel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SegPanel (see VARARGIN)

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

% Choose default command line output for SegPanel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes projSetupGUI wait for user response (see UIRESUME)
%uiwait(handles.output);

initialize_gui(hObject, handles, false);

% UIWAIT makes SegPanel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SegPanel_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function density_Callback(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density as text
%        str2double(get(hObject,'String')) returns contents of density as a double
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.metricdata.density = density;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function volume_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function volume_Callback(hObject, eventdata, handles)
% hObject    handle to volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of volume as text
%        str2double(get(hObject,'String')) returns contents of volume as a double
volume = str2double(get(hObject, 'String'));
if isnan(volume)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new volume value
handles.metricdata.volume = volume;
guidata(hObject,handles)

% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mass = handles.metricdata.density * handles.metricdata.volume;
set(handles.mass, 'String', mass);

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

initialize_gui(gcbf, handles, true);

% --------------------------------------------------------------------
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to unitgroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (hObject == handles.english)
    set(handles.text4, 'String', 'lb/cu.in');
    set(handles.text5, 'String', 'cu.in');
    set(handles.text6, 'String', 'lb');
else
    set(handles.text4, 'String', 'kg/cu.m');
    set(handles.text5, 'String', 'cu.m');
    set(handles.text6, 'String', 'kg');
end

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, segpanel means
% we are we are just re-initializing a GUI by calling segpanel from the cmd line
% while segpanel is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

handles.metricdata.density = 0;
handles.metricdata.volume  = 0;

%set(handles.density, 'String', handles.metricdata.density);
% set(handles.volume,  'String', handles.metricdata.volume);
% set(handles.mass, 'String', 0);
% 
% set(handles.unitgroup, 'SelectedObject', handles.english);
% 
% set(handles.text4, 'String', 'lb/cu.in');
% set(handles.text5, 'String', 'cu.in');
% set(handles.text6, 'String', 'lb');

% Update handles structure
guidata(handles.figure1, handles);





% --- Executes on button press in pushLoadStacks.
function pushLoadStacks_Callback(hObject, eventdata, handles)
% hObject    handle to pushLoadStacks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%allow user to choose directory
try
    basePath = uigetdir('Please select directory of movies to be collected');
    if basePath == 0
        disp('Interrupted by user.');
        return;
    end
end

%find all makiAnalysis files in chosen directory
fileList = searchFiles('tif',[],basePath,1);
channelDir = unique(fileList(:,2));
if strcmp(channelDir{1},pwd) == 1
   channelDir(1) = [];
end
    
%allow user to choose files
selectIdx = listSelectGUI(channelDir,[],'move');

%get number of chosen files
channelDirs = channelDir(selectIdx,:);
numChannels = length(selectIdx);

if numChannels == 0 %no directory has been specified
    f = warndlg('Please specify at least one channel', 'Warning', 'modal');
    return;
end

clear fileList
for n_idx = 1:numChannels
    fileList{n_idx} = searchFiles('tif',[],channelDirs{n_idx},0);    
    numImages(n_idx) = length(fileList{n_idx});
end

if numChannels > 1
    if  length(unique(numImages)) > 1
        f = warndlg(['The selected directores have different number of images.  Please check again before proceed.'], 'Warning', 'modal');
        return;
    end
    for c_idx = 1:numChannels
        [this_row this_col] = size(imread([fileList{c_idx}{1,2} filesep fileList{c_idx}{1,1}]));
        r(c_idx) = this_row;
        c(c_idx) = this_col;
    end
    if sum(abs(diff(r))) > 0 | sum(abs(diff(c))) > 0
        f = warndlg(['The selected directores have images with different dimensions.  Please check again before proceed.'], 'Warning', 'modal');
        return;        
    end
    
    slider_channel_step(1) = 1/(numChannels - 1);
    slider_channel_step(2) = 1/(numChannels - 1);
    set(handles.slider_channel, 'Visible', 'on','Max',numChannels,...
        'Value',1,'Min',1,'sliderstep',slider_channel_step);    
else
    set(handles.slider_channel, 'Visible', 'off');
end
    
mask_idx = 1;
slider_time_step(1) = 1/numImages(1);
slider_time_step(2) = 1/numImages(1);

channel_idx = 1;
thismask = imread([fileList{channel_idx}{mask_idx,2} filesep fileList{channel_idx}{mask_idx,1}]);
set(handles.slider_time, 'Visible', 'on','Max',numImages(1),...
    'Value',1,'Min',1,'sliderstep',slider_time_step);
 
axes(handles.axes_image);    imshow(thismask, []);   colormap(handles.colormap_type);
handles.thismask = thismask;
handles.mask_idx = mask_idx;

handles.rowSize = size(thismask,1);
handles.colSize = size(thismask,2);
handles.numChannels = numChannels;
handles.numImages = numImages;
handles.channel_idx = channel_idx;
handles.fileList = fileList;
handles.channelDirs = channelDirs;

%text(handles.colSize+5,20,['Image Index:' num2str(mask_idx)],'BackgroundColor',[.7 .9 .7]);
set(handles.text_raw, 'String', ['Region Of Interest: Image ' num2str(mask_idx) ' in Channel ' num2str(channel_idx)]);    

set(handles.text_raw , 'Visible', 'on');
set(handles.push_applySeg , 'Visible', 'on');
set(handles.panel_seg, 'Visible', 'on');
set(handles.push_segagain, 'Visible', 'off');

Closure_Size_Max = 100;
slider_channel_step(1) = 1/(Closure_Size_Max-1);
slider_channel_step(2) = 1/(Closure_Size_Max-1);
set(handles.slider_closure, 'Visible', 'on','Max',Closure_Size_Max,...
    'Value',5,'Min',0,'sliderstep',slider_channel_step);
handles.closure_size = 5;
handles.segType = 4;
handles.isseg = 0;
% Update handles structure
guidata(hObject, handles);

update_Seg(handles);

% --- Executes on button press in backbutton.
function backbutton_Callback(hObject, eventdata, handles)
% hObject    handle to backbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%delete(hObject);
delete(handles.figure1);
%uiresume(handles.output);




function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_applySeg.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to push_applySeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in push_play.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to push_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
                                                                                                                                                              

% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in push_crop2.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to push_crop2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in push_crop3.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to push_crop3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in push_crop4.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to push_crop4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton2


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','off');
% Hint: place code in OpeningFcn to populate axes4




% --- Executes during object creation, after setting all properties.
function axes_segimage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_segimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','off');
% Hint: place code in OpeningFcn to populate axes_segimage


                                                                                                                                                        

% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','off');
% Hint: place code in OpeningFcn to populate axes5


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



% --- Executes on slider movement.
function slider_time_Callback(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.mask_idx = get(handles.slider_time,'Value');
handles.mask_idx = round(handles.mask_idx);
set(handles.slider_time,'Value',handles.mask_idx);
thismask = imread([handles.fileList{handles.channel_idx}{handles.mask_idx,2} filesep handles.fileList{handles.channel_idx}{handles.mask_idx,1}]);
axes(handles.axes_image);
imshow(thismask, []);   colormap(handles.colormap_type);
set(handles.text_raw, 'String', ['Region Of Interest: Image ' num2str(handles.mask_idx) ' in Channel ' num2str(handles.channel_idx)]);    

handles.thismask = thismask;

if handles.isseg == 1
    thissegmask = imread(handles.segfileList{handles.channel_idx}{handles.mask_idx});
    axes(handles.axes_segimage);
    imshow(thissegmask, []);   colormap(handles.colormap_type);
    set(handles.text_raw, 'String', ['Segmentation: Image ' num2str(handles.mask_idx) ' in Channel ' num2str(handles.channel_idx)]);
else
    update_Seg(handles);
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on slider movement.
function slider_channel_Callback(hObject, eventdata, handles)
% hObject    handle to slider_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.channel_idx = get(handles.slider_channel,'Value');
handles.channel_idx = round(handles.channel_idx);
set(handles.slider_channel, 'String', num2str(handles.channel_idx));

thismask = imread([handles.fileList{handles.channel_idx}{handles.mask_idx,2} filesep handles.fileList{handles.channel_idx}{handles.mask_idx,1}]);
axes(handles.axes_image);
imshow(thismask, []);   colormap(handles.colormap_type);
set(handles.text_raw, 'String', ['Region Of Interest: Image ' num2str(handles.mask_idx) ' in Channel ' num2str(handles.channel_idx)]);    
handles.thismask = thismask;

if handles.isseg == 1
    thissegmask = imread(handles.segfileList{handles.channel_idx}{handles.mask_idx});
    axes(handles.axes_segimage);
    imshow(thissegmask, []);   colormap(handles.colormap_type);
    set(handles.text_raw, 'String', ['Segmentation: Image ' num2str(handles.mask_idx) ' in Channel ' num2str(handles.channel_idx)]);
else
    update_Seg(handles);
end

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in push_applySeg.
function push_applySeg_Callback(hObject, eventdata, handles)
% hObject    handle to push_applySeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

thispath=uigetdir('','Select an output directory');
if thispath==0
    guidata(hObject, handles);    
    thispath = warndlg('Please specify an output directory then press any crop button again', 'Warning', 'modal');
    return;
end

for c_idx=1:handles.numChannels
    [fpath,dirName,fno,fext]=getFilenameBody(handles.channelDirs{c_idx});
    if ~exist([thispath filesep dirName num2str(fno)])    
        mkdir(thispath,[dirName num2str(fno)]);
    end
    segDir{c_idx} = [dirName num2str(fno)];
end
cellmaskDir = [thispath filesep 'cell_mask'];
if ~exist(cellmaskDir)
    mkdir(cellmaskDir);
end

edgecellDir = [thispath filesep 'edge_cell'];
if ~exist(edgecellDir)
    mkdir(edgecellDir);
end

% ROImaskDir = [thispath filesep 'ROI-bgsub'];
% 
% if ~exist(ROImaskDir)
%     mkdir(ROImaskDir);
% end

% Initializing waitbar
h=waitbar(0,'Processing...');

% Processing files
for n_idx=1:handles.numImages(1)
    thismask = imread([handles.fileList{handles.channel_idx}{n_idx,2} filesep handles.fileList{handles.channel_idx}{n_idx,1}]);
    switch handles.segType
        case 1
            % Code for when radiobutton_seg1 is selected.
            % Otsu Thresholding
            segmask = otsuSeg(thismask, handles.closure_size);
        case 2
            % Code for when radiobutton_seg2 is selected.
            % Median Thresholding
            segmask = medianSeg(thismask, handles.closure_size);
        case 3
            % Code for when radiobutton_seg3 is selected.
            % Phase Contrast Segmentation
            segmask = phasecontrastSeg(thismask, handles.closure_size);
        case 4
            % Code for when radiobutton_seg4 is selected.
            % Customized Segmentation
            segmask = CustomizedSeg(thismask, handles.closure_size);            
        otherwise
            % Code for when there is no match.
    end
    filename = [cellmaskDir filesep 'mask_' handles.fileList{handles.channel_idx}{n_idx,1}];
    % Write file to disk
    imwrite(segmask,filename, 'Compression', 'none');

%    ROImaskfilename = [ROImaskDir filesep handles.fileList{handles.channel_idx}{n_idx,1}];
%    imwrite(uint16(segmask) .* uint16(thismask), ROImaskfilename, 'Compression', 'none');    

    for c_idx=1:handles.numChannels
        % Prepare filename with path
        thismask = imread([handles.fileList{c_idx}{n_idx,2} filesep handles.fileList{c_idx}{n_idx,1}]);
        filename=[thispath,filesep,segDir{c_idx},filesep,handles.fileList{c_idx}{n_idx,1},'.tif'];
        handles.segfileList{c_idx}{n_idx} = filename;

%         c = contourc(double(segmask), [0 0]);
%         pixel_list = round([c(1,2:end)' c(2,2:end)']);
        myedge = bwperim(segmask);
        [my_x, my_y] = find(myedge == 1);
        pixel_list = [my_x, my_y];
        h_prot_control = figure('Visible','Off');
        im3color = zeros(size(thismask));
        im3color(:,:,1) = mat2gray(thismask);
        im3color(:,:,2) = mat2gray(thismask);
        im3color(:,:,3) = mat2gray(thismask);        
        for i=1:size(pixel_list,1)
            im3color(pixel_list(i,1),pixel_list(i,2),1) = 1;
        end
        imwrite(im3color,filename, 'Compression', 'none');        
    end
    
    for c_idx=1:1
        % Prepare filename with path
        thismask = imread([handles.fileList{c_idx}{n_idx,2} filesep handles.fileList{c_idx}{n_idx,1}]);
        filename=[edgecellDir,filesep,'img_edge_',handles.fileList{c_idx}{n_idx,1},'.tif'];
        handles.segfileList{c_idx}{n_idx} = filename;
        myedge = bwperim(segmask);
        [my_x, my_y] = find(myedge == 1);
        pixel_list = [my_x, my_y];
        h_prot_control = figure('Visible','Off');
        im3color = zeros(size(thismask));
        im3color(:,:,1) = mat2gray(thismask);
        im3color(:,:,2) = mat2gray(thismask);
        im3color(:,:,3) = mat2gray(thismask);        
        for i=1:size(pixel_list,1)
            im3color(pixel_list(i,1),pixel_list(i,2),1) = 1;
        end
        imwrite(im3color,filename, 'Compression', 'none');        
    end    
    % Update waitbar
    waitbar(n_idx/handles.numImages(1),h);
end
close(h);
handles.isseg = 1;
handles.thispath =thispath;
set(handles.push_play, 'Visible', 'on');
set(handles.push_segagain, 'Visible', 'on');
set(handles.panel_seg, 'Visible', 'off');
set(handles.pushLoadStacks, 'Visible', 'off');
set(handles.push_applySeg, 'Visible', 'off');
set(handles.bgRemoval_button , 'Visible', 'on');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in push_play.
function push_play_Callback(hObject, eventdata, handles)
% hObject    handle to push_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

channel_idx = handles.channel_idx;
for mask_idx = 1:handles.numImages

    thismask = imread([handles.fileList{channel_idx}{mask_idx,2} filesep handles.fileList{channel_idx}{mask_idx,1}]);
    axes(handles.axes_image);    imshow(thismask, []);   colormap(handles.colormap_type);

    thiscropmask = imread([handles.segfileList{channel_idx}{mask_idx}]);
    axes(handles.axes_segimage);    imshow(thiscropmask, []);   colormap(handles.colormap_type);

    set(handles.slider_time,'Value',mask_idx);

    set(handles.text_cropped, 'String', ['Segmentation: Image ' num2str(mask_idx) ' in Channel ' num2str(channel_idx)]);
    set(handles.text_raw, 'String', ['Region Of Interest: Image ' num2str(mask_idx) ' in Channel ' num2str(channel_idx)]);
    pause(0.02);
end
handles.mask_idx = mask_idx;

guidata(hObject, handles);

% --- Executes on button press in push_crop2.
function push_crop2_Callback(hObject, eventdata, handles)
% hObject    handle to push_crop2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in push_crop3.
function push_crop3_Callback(hObject, eventdata, handles)
% hObject    handle to push_crop3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in push_crop4.
function push_crop4_Callback(hObject, eventdata, handles)
% hObject    handle to push_crop4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function axes_logo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_logo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_logo

imagePath=['penguin.bmp'];
axes(hObject);
imshow(imread(imagePath), []);                                                                                                                          



% --- Executes on button press in push_segagain.
function push_segagain_Callback(hObject, eventdata, handles)
% hObject    handle to push_segagain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.push_segagain, 'Visible', 'off');
set(handles.panel_seg, 'Visible', 'on');
set(handles.push_play, 'Visible', 'off');
set(handles.pushLoadStacks, 'Visible', 'on');
set(handles.push_applySeg, 'Visible', 'on');
set(handles.bgRemoval_button , 'Visible', 'off');

handles.isseg = 0;
guidata(hObject, handles);


% --- Executes when selected object is changed in panel_seg.
function panel_seg_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_seg 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject,'Tag')   % Get Tag of selected object
    case 'radiobutton_seg1'
        % Code for when radiobutton_seg1 is selected.
        % Otsu Thresholding
        handles.segType = 1;
    case 'radiobutton_seg2'
        % Code for when radiobutton_seg2 is selected.
        % Median Thresholding
        handles.segType = 2;        
    case 'radiobutton_seg3'
        % Code for when radiobutton_seg3 is selected.
        % Phase Contrast Segmentation
        handles.segType = 3;        
    case 'radiobutton_seg4'
        % Code for when radiobutton_seg4 is selected.
        % Customized Segmentation
        handles.segType = 4;                
    otherwise
        % Code for when there is no match.
        handles.segType = 1;
end

update_Seg(handles);

% Update handles structure
guidata(hObject, handles);


function update_Seg(handles)

thismask = double(handles.thismask);

switch handles.segType   % Get Tag of selected object
    case 1
        % Code for when radiobutton_seg1 is selected.
        % Otsu Thresholding
        segmask = otsuSeg(thismask, handles.closure_size);
    case 2
        % Code for when radiobutton_seg2 is selected.
        % Median Thresholding
        segmask = medianSeg(thismask, handles.closure_size);
    case 3
        % Code for when radiobutton_seg3 is selected.
        % Phase Contrast Segmentation
        segmask = phasecontrastSeg(thismask, handles.closure_size);
    case 4
        % Code for when radiobutton_seg4 is selected.
        % Customized Segmentation
        segmask = CustomizedSeg(thismask, handles.closure_size);        
    otherwise
        % Code for when there is no match.
end

axes(handles.axes_segimage);
%imshow(segmask, []);
c = contourc(double(segmask), [0 0]);
pixel_list = [c(1,2:end)' c(2,2:end)'];
imshow(handles.thismask, []);
colormap(handles.colormap_type);
hold on; plot(pixel_list(:,1), pixel_list(:,2), 'r.','MarkerSize',2); hold off;



% --- Executes on slider movement.
function slider_closure_Callback(hObject, eventdata, handles)
% hObject    handle to slider_closure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.closure_size = get(handles.slider_closure,'Value');
handles.closure_size = round(handles.closure_size);
set(handles.slider_closure,'Value',handles.closure_size);
set(handles.text_closure_size,'String',handles.closure_size);

update_Seg(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider_closure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_closure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
manualFreeHand;



% --- Executes on button press in bgRemoval_button.
function bgRemoval_button_Callback(hObject, eventdata, handles)
% hObject    handle to bgRemoval_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ROImaskDir = [handles.thispath filesep 'ROI-bgsub'];
maskDir = [handles.thispath filesep 'masks'];

if ~exist(ROImaskDir)
    mkdir(ROImaskDir);
end

% Initializing waitbar
h=waitbar(0,'Processing...');

% Processing files
for n_idx=1:handles.numImages(1)
    maskfilename = [maskDir filesep handles.fileList{handles.channel_idx}{n_idx,1}];
%    maskfilename = [maskDir filesep handles.fileList{handles.channel_idx}{1,1}];
    segmask = imread(maskfilename);
    ROImaskfilename = [ROImaskDir filesep handles.fileList{handles.channel_idx}{n_idx,1}];

    for c_idx=1:handles.numChannels
        % Prepare filename with path
        thismask = imread([handles.fileList{c_idx}{n_idx,2} filesep handles.fileList{c_idx}{n_idx,1}]);
        imwrite(uint16(segmask).*uint16(thismask), ROImaskfilename, 'Compression', 'none');    
    end
    % Update waitbar
    waitbar(n_idx/handles.numImages(1),h);
end
close(h);

