function varargout = SignalPlotting(varargin)
% SIGNALPLOTTING MATLAB code for SignalPlotting.fig
%      SIGNALPLOTTING, by itself, creates a new SIGNALPLOTTING or raises the existing
%      singleton*.
%
%      H = SIGNALPLOTTING returns the handle to a new SIGNALPLOTTING or the handle to
%      the existing singleton*.
%
%      SIGNALPLOTTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIGNALPLOTTING.M with the given input arguments.
%
%      SIGNALPLOTTING('Property','Value',...) creates a new SIGNALPLOTTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SignalPlotting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SignalPlotting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% edit_celltrackmatfilename the above text to modify the response to help SignalPlotting

% Last Modified by GUIDE v2.5 07-Apr-2013 14:34:09

% Begin initialization code - DO NOT EDIT_CELLTRACKMATFILENAME
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SignalPlotting_OpeningFcn, ...
                   'gui_OutputFcn',  @SignalPlotting_OutputFcn, ...
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
% End initialization code - DO NOT EDIT_CELLTRACKMATFILENAME


% --- Executes just before SignalPlotting is made visible.
function SignalPlotting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SignalPlotting (see VARARGIN)

% Choose default command line output for SignalPlotting
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SignalPlotting wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SignalPlotting_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_signalmatfilename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_signalmatfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_signalmatfilename as text
%        str2double(get(hObject,'String')) returns contents of edit_signalmatfilename as a double


% --- Executes during object creation, after setting all properties.
function edit_signalmatfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_signalmatfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit_celltrackmatfilename controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_choosefile.
function pushbutton_choosefile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_choosefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname,FilterIndex] = uigetfile(...
{'*.mat','cell track mat file(*.mat)';
   '*.*',  'All Files (*.*)'}, 'Pick a file');
if FilterIndex~=0
    set(handles.edit_signalmatfilename,'String',filename);
    handles.signalfilename = fullfile(pathname, filename);
end
guidata(hObject, handles);  



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname,FilterIndex] = uigetfile(...
{'*.mat','cell track mat file(*.mat)';
   '*.*',  'All Files (*.*)'}, 'Pick a file');
if FilterIndex~=0
    set(handles.edit_celltrackmatfilename,'String',filename);
    handles.signalfilename = fullfile(pathname, filename);
end
guidata(hObject, handles);  

% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load(get(handles.edit_signalmatfilename,'String'));
s = whos;
cind = 1;
for i=1:length(s)
    varsize = s(i).size;
    if varsize(1)~=0 & varsize(1)~=1
        newString{cind} = s(i).name;
        eval(['handles.' s(i).name '=' s(i).name ';']);
        cind=cind+1;
        clear(s(i).name);
    end
end
%handles.timestep = timestep;
%clear timestep;
set(handles.popupmenu_signal,'String',newString);
handles.selected_cells = selected_cells;
load(get(handles.edit_celltrackmatfilename,'String'));
handles.cellpath = cellpath;
handles.sisterList = sisterList;
handles.bg = bg;
handles.mark = mark;
guidata(hObject, handles); 




% --- Executes on selection change in popupmenu_stagePos.
function popupmenu_stagePos_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_stagePos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.edit_stageInfo,'String',handles.stageName{get(hObject,'Value')});
tokens   = regexp(handles.stageName{get(hObject,'Value')}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
if ~isempty(tokens)
    row = tokens{1}{1};
    col = tokens{1}{2};
    set(handles.edit_row,'String',row);
    set(handles.edit_col,'String',col);
else
    set(handles.edit_row,'String',num2str(get(hObject,'Value')));
    set(handles.edit_col,'String','1');
    row = get(hObject,'Value');
    col = 1;
end

field = get(handles.edit_field,'String');

H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
signal_name = ['/field' num2str(field) '/signals'];
timestamp_name = ['/field' num2str(field) '/timestamp'];
selectedcells_name = ['/field' num2str(field) '/selectedcells'];


if exist(fullfile(handles.SourceF,H5filename),'file')
    fid = H5F.open(fullfile(handles.SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if H5L.exists(fid,signal_name,'H5P_DEFAULT')
        H5F.close(fid);
        signalinfo = h5info(fullfile(handles.SourceF,H5filename), signal_name);
        startind = double([1 1 1]);
        countind = [signalinfo.Dataspace.Size(1) signalinfo.Dataspace.Size(2) 4];
        signals = h5read(fullfile(handles.SourceF,H5filename),signal_name,startind, countind);
        timestamp = h5read(fullfile(handles.SourceF,H5filename),timestamp_name);
        handles.signals = signals; %  CellNo,Time,Signal
        handles.timestamp = timestamp;
        
    end
    fid = H5F.open(fullfile(handles.SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if H5L.exists(fid,selectedcells_name,'H5P_DEFAULT')
        H5F.close(fid);
        
        selected_cells = h5read(fullfile(handles.SourceF,H5filename),selectedcells_name);
        handles.selected_cells=selected_cells;
    end
    
    if ~isempty(selected_cells)
        for scell = selected_cells'
            PosInd = find(signals(scell,:,1));
            axes(handles.axes1);
            if scell==selected_cells(1)
                plot(timestamp(PosInd)/60,signals(scell,PosInd,1),'color',[0.7 0.7 0.7]);
            else
                hold on;
                plot(timestamp(PosInd)/60,signals(scell,PosInd,1),'color',[0.7 0.7 0.7]);
                hold off;
                
            end
        end
        
        for tp=1:length(timestamp)
            PosInd = find(signals(:,tp,1));
            median_Signal(tp) = median(signals(PosInd,tp,1));
        end
        axes(handles.axes1);
        hold on;
        plot(timestamp/60,median_Signal,'k');
        hold off;
    end
    
    %signalNames = get(handles.popupmenu_regionVar1,'String');
    
else
    set(handles.edit_commu,'String','Check to make sure that H5 file exists');
    return;
end

% cind = 1;
% for i=1:length(s)
%     varsize = s(i).size;
%     if varsize(1)~=0 & varsize(1)~=1
%         newString{cind} = s(i).name;
%         cind=cind+1;
%         clear(s(i).name);
%     end
% end
% 
% set(handles.popupmenu_signal,'String',newString);
guidata(hObject, handles);  


% --- Executes on selection change in popupmenu_signal.
function popupmenu_signal_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_signal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_signal


scell = str2num(get(handles.edit_cellNo,'String'));
PosTime = find(signal(:,scell,1));
figure;
signalNames = get(handles.popupmenu_regionVar1,'String');
hold on;plot(timestamp(PosTime)/60,signal(PosTime,scell,1)/median(signal(PosTime,scell,1)),'b');
s1Loc = get(handles.popupmenu_regionVar1,'Value');
legendList{1} = signalNames{s1Loc};
xlabel('Time(hour)');


% --- Executes on button press in pushbutton_plot.
function pushbutton_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));



H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
signal_name = ['/field' num2str(field) '/' get(handles.edit_outputname,'String')];
timestamp_name = ['/field' num2str(field) '/timestamp'];



if exist(fullfile(handles.SourceF,H5filename),'file')
    fid = H5F.open(fullfile(handles.SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if H5L.exists(fid,signal_name,'H5P_DEFAULT')
        H5F.close(fid);
        legendList=cell(1);
        signalinfo = h5info(fullfile(handles.SourceF,H5filename), signal_name);
        startind = double([1 1 1]);
        countind = [signalinfo.Dataspace.Size(1) signalinfo.Dataspace.Size(2) 4];
        signal = permute(h5read(fullfile(handles.SourceF,H5filename),signal_name,startind, countind),[2 1 3]);
        timestamp = double(h5read(fullfile(handles.SourceF,H5filename),timestamp_name));
        
        scell = str2num(get(handles.edit_cellNo,'String'));
        PosTime = find(signal(:,scell,1));
        figure;
        signalNames = get(handles.popupmenu_regionVar1,'String');
        hold on;plot(timestamp(PosTime)/60,signal(PosTime,scell,1)/median(signal(PosTime,scell,1)),'b');
        s1Loc = get(handles.popupmenu_regionVar1,'Value');
        legendList{1} = signalNames{s1Loc};
        if get(handles.checkbox_variable2,'Value')
            hold on;plot(timestamp(PosTime)/60,signal(PosTime,scell,2)/median(signal(PosTime,scell,2)),'r');
            s2Loc = get(handles.popupmenu_regionVar2,'Value');
            legendList{2} = signalNames{s2Loc};
        end
        if get(handles.checkbox_variable3,'Value')
            hold on;plot(timestamp(PosTime)/60,signal(PosTime,scell,3)/median(signal(PosTime,scell,3)),'g');
            s3Loc = get(handles.popupmenu_regionVar3,'Value');
            legendList{3} =  signalNames{s3Loc};
        end
        if get(handles.checkbox_variable4,'Value')
            hold on;plot(timestamp(PosTime)/60,signal(PosTime,scell,4)/median(signal(PosTime,scell,4)),'k');
            s4Loc = get(handles.popupmenu_regionVar4,'Value');
            legendList{4} = signalNames{s4Loc};
        end
        legend(legendList);
        xlabel('Time(hour)');
    else
        set(handles.edit_commu,'String',[signal_name ' does not exist']);
        return;
        
    end

else
    set(handles.edit_commu,'String','Check to make sure that H5 file exists');
    return;
end

% --- Executes during object creation, after setting all properties.
function popupmenu_signal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_datatype.
function popupmenu_datatype_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_datatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_datatype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_datatype


% --- Executes during object creation, after setting all properties.
function popupmenu_datatype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_datatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plottype.
function popupmenu_plottype_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plottype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plottype


% --- Executes during object creation, after setting all properties.
function popupmenu_plottype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_framerange_Callback(hObject, eventdata, handles)
% hObject    handle to edit_framerange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_framerange as text
%        str2double(get(hObject,'String')) returns contents of edit_framerange as a double


% --- Executes during object creation, after setting all properties.
function edit_framerange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_framerange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit_celltrackmatfilename controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_normalization.
function popupmenu_normalization_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_normalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_normalization contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_normalization


% --- Executes during object creation, after setting all properties.
function popupmenu_normalization_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_normalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_celltrackmatfilename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_celltrackmatfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_celltrackmatfilename as text
%        str2double(get(hObject,'String')) returns contents of edit_celltrackmatfilename as a double


% --- Executes during object creation, after setting all properties.
function edit_celltrackmatfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_celltrackmatfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit_celltrackmatfilename controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minI as text
%        str2double(get(hObject,'String')) returns contents of edit_minI as a double


% --- Executes during object creation, after setting all properties.
function edit_minI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxI as text
%        str2double(get(hObject,'String')) returns contents of edit_maxI as a double


% --- Executes during object creation, after setting all properties.
function edit_maxI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_commonNorm.
function checkbox_commonNorm_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_commonNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_commonNorm



function edit_timeframe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_timeframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_timeframe as text
%        str2double(get(hObject,'String')) returns contents of edit_timeframe as a double


% --- Executes during object creation, after setting all properties.
function edit_timeframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_timeframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_movingwindow_Callback(hObject, eventdata, handles)
% hObject    handle to edit_movingwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_movingwindow as text
%        str2double(get(hObject,'String')) returns contents of edit_movingwindow as a double


% --- Executes during object creation, after setting all properties.
function edit_movingwindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_movingwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_smoothingLogic.
function checkbox_smoothingLogic_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_smoothingLogic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_smoothingLogic


% --- Executes on selection change in popupmenu_colormap.
function popupmenu_colormap_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_colormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_colormap


% --- Executes during object creation, after setting all properties.
function popupmenu_colormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_rank.
function checkbox_rank_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_rank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_rank



function edit_clusterno_Callback(hObject, eventdata, handles)
% hObject    handle to edit_clusterno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_clusterno as text
%        str2double(get(hObject,'String')) returns contents of edit_clusterno as a double


% --- Executes during object creation, after setting all properties.
function edit_clusterno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_clusterno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_frameduration_Callback(hObject, eventdata, handles)
% hObject    handle to edit_frameduration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_frameduration as text
%        str2double(get(hObject,'String')) returns contents of edit_frameduration as a double


% --- Executes during object creation, after setting all properties.
function edit_frameduration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_frameduration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_refframe2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_refframe2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_refframe2 as text
%        str2double(get(hObject,'String')) returns contents of edit_refframe2 as a double


% --- Executes during object creation, after setting all properties.
function edit_refframe2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_refframe2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_refframe1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_refframe1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_refframe1 as text
%        str2double(get(hObject,'String')) returns contents of edit_refframe1 as a double


% --- Executes during object creation, after setting all properties.
function edit_refframe1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_refframe1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxnormI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxnormI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxnormI as text
%        str2double(get(hObject,'String')) returns contents of edit_maxnormI as a double


% --- Executes during object creation, after setting all properties.
function edit_maxnormI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxnormI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_row_Callback(hObject, eventdata, handles)
% hObject    handle to edit_row (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_row as text
%        str2double(get(hObject,'String')) returns contents of edit_row as a double


% --- Executes during object creation, after setting all properties.
function edit_row_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_row (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_col_Callback(hObject, eventdata, handles)
% hObject    handle to edit_col (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_col as text
%        str2double(get(hObject,'String')) returns contents of edit_col as a double


% --- Executes during object creation, after setting all properties.
function edit_col_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_col (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_field_Callback(hObject, eventdata, handles)
% hObject    handle to edit_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_field as text
%        str2double(get(hObject,'String')) returns contents of edit_field as a double


% --- Executes during object creation, after setting all properties.
function edit_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_plane_Callback(hObject, eventdata, handles)
% hObject    handle to edit_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_plane as text
%        str2double(get(hObject,'String')) returns contents of edit_plane as a double


% --- Executes during object creation, after setting all properties.
function edit_plane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ndfilename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ndfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ndfilename as text
%        str2double(get(hObject,'String')) returns contents of edit_ndfilename as a double
handles.ndfilename = get(hObject,'String');
guidata(hObject, handles);  

% --- Executes during object creation, after setting all properties.
function edit_ndfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ndfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_locatendfile.
function pushbutton_locatendfile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_locatendfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,SourceF,FilterIndex] = uigetfile('*.nd', 'Choose metamorph ND file','C:\computation\02-03-2013\02032013-r1.nd');
if FilterIndex~=0
    set(handles.edit_ndfilename,'String',filename);
    handles.ndfilename = filename;
    
    handles.SourceF = SourceF;
    set(handles.edit_sourceF,'String',SourceF);
end

guidata(hObject, handles);  


% --- Executes on button press in pushbutton_loadbyndfile.
function pushbutton_loadbyndfile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadbyndfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[notp stagePos stageName waveName] = readndfile(fullfile(handles.SourceF,handles.ndfilename));

if notp==-1
    return;
end

set(handles.popupmenu_stagePos,'String',stagePos);
set(handles.popupmenu_stagePos,'Value',1);
set(handles.edit_stageInfo,'String',stageName{1});
handles.stageName = stageName;
handles.channelnames = waveName;
handles.prefix = handles.ndfilename(1:(end-3));
tokens   = regexp(stageName{1}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
if ~isempty(tokens)
    
    row = tokens{1}{1};
    col = tokens{1}{2};
    set(handles.edit_row,'String',row);
    set(handles.edit_col,'String',col);
else
    set(handles.edit_row,'String','1');
    set(handles.edit_col,'String','1');
end

guidata(hObject, handles);  



function [notp stagePos stageName waveName] = readndfile(filename)
% Search for number of string matches per line.  
notp=-1;
stagePos = [];
stageName = [];
waveName = [];


if exist(filename,'file')
    fid = fopen(filename);
    y = 0;
    tline = fgetl(fid);
    sind = 1;
    wind = 1;
    notp=0;
    while ischar(tline)
        
        % Find number of time points
        
        testInd = regexp(tline,'NTimePoints');
        num = length(testInd);
        if num > 0
            tp  = regexp(tline, '(?<="NTimePoints", )\d+', 'match');
            notp = str2num(tp{1});
        end
        
        
        % Find stage naming
        testInd = regexp(tline,'Stage\d+');
        num = length(testInd);
        if num > 0
            stage  = regexp(tline, '(?<=")\w+(?=",)', 'match');
            stagePos{sind,1} = stage{1};
            stagename  = regexp(tline, '(?<="Stage\d+", ")\w+(?=")', 'match');
            stageName{sind,1} = stagename{1};
            sind=sind+1;
        end
        
        % Find stage naming
        testInd = regexp(tline,'WaveName\d+');
        num = length(testInd);
        if num > 0
            wavename1  = regexp(tline, '(?<="WaveName\d+", ")\w+(?=_)', 'match');
            wavename2  = regexp(tline, '(?<="WaveName\d+", "\w+_)\w+(?=")', 'match');
            waveName{wind} = ['w' num2str(wind) wavename1{1} '-' wavename2{1}];
            wind=wind+1;
        end
        
        tline = fgetl(fid);
    end
    fclose(fid);
end

function edit_sourceF_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sourceF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sourceF as text
%        str2double(get(hObject,'String')) returns contents of edit_sourceF as a double
handles.SourceF = get(hObject,'String');
guidata(hObject, handles);  

% --- Executes during object creation, after setting all properties.
function edit_sourceF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sourceF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function popupmenu_stagePos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_stagePos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stageInfo_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stageInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stageInfo as text
%        str2double(get(hObject,'String')) returns contents of edit_stageInfo as a double


% --- Executes during object creation, after setting all properties.
function edit_stageInfo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stageInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_commu_Callback(hObject, eventdata, handles)
% hObject    handle to edit_commu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_commu as text
%        str2double(get(hObject,'String')) returns contents of edit_commu as a double


% --- Executes during object creation, after setting all properties.
function edit_commu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_commu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
