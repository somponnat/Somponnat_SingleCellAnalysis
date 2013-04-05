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

% Last Modified by GUIDE v2.5 24-Jan-2013 17:02:19

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

% --- Executes on selection change in popupmenu_signal.
function popupmenu_signal_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_signal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_signal

% --- Executes on button press in pushbutton_plot.
function pushbutton_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sisterList = handles.sisterList;
cellpath = handles.cellpath;
varind = get(handles.popupmenu_signal,'Value');
varnames = get(handles.popupmenu_signal,'String');
varname = varnames{varind};
eval(['myvar=' 'handles.' varname ';']);

timeframes = str2num(get(handles.edit_timeframe,'String'));

switch get(handles.popupmenu_datatype,'Value')
    case 1 %all
        plotind = handles.selected_cells;
        pairdata = 0;
    case 2 %non-dividing
        plotind = find(sisterList{end}(:,1)==-1);
        pairdata = 0;
    case 3 %dividing
        sind = find(sisterList{end}(:,1)~=-1);
        cind = 1;
        plotind = [];
        for sl=1:length(sind)
            tempind = find(handles.selected_cells == sind(sl));
            if isempty(tempind)
                plotind(cind) = sind(sl);
                plotind(cind+1) = sisterList{end}(sind(sl),1);
                cind = cind+2;
            end
        end
        pairdata = 1;
        
end
plotInd=1;
switch get(handles.popupmenu_normalization,'Value')
    case 1
        if pairdata==1
            for c=plotind
                if get(handles.checkbox_commonNorm,'Value')
                    sigmax = str2num(get(handles.edit_maxI,'String'));
                    sigmin = str2num(get(handles.edit_minI,'String'));
                else
                    sigmax = max(myvar{c}(timeframes,2));
                    posInds = find(myvar{c}(timeframes,2)>0);
                    sigmin = min(myvar{c}(posInds,2));
                end
                normvar(:,c) = (myvar{c}(timeframes,2)-sigmin)/(sigmax-sigmin);
                normvar(:,c+1) = (myvar{c+1}(timeframes,2)-sigmin)/(sigmax-sigmin);
            end
        else
            for c=plotind
                if get(handles.checkbox_commonNorm,'Value')
                    sigmax = str2num(get(handles.edit_maxI,'String'));
                    sigmin = str2num(get(handles.edit_minI,'String'));
                    normvar(:,plotInd) = (myvar(timeframes,c)-sigmin)/(sigmax-sigmin);
                    plotInd=plotInd+1;
                else
                    sigmax = max(myvar{c}(timeframes,2));
                    posInds = find(myvar{c}(timeframes,2)>0);
                    sigmin = min(myvar{c}(posInds,2));
                    normvar(:,plotInd) = (myvar{c}(timeframes,2)-sigmin)/(sigmax-sigmin);
                    negInds = find(normvar(:,plotInd)<0);
                    normvar(negInds,plotInd) = 0;
                    plotInd=plotInd+1;
                end
            end
        end
    case 2
        if pairdata==1
            for c=plotind
                if get(handles.checkbox_commonNorm,'Value')
                    sigmax = str2num(get(handles.edit_maxI,'String'));
                    sigmin = str2num(get(handles.edit_minI,'String'));
                else
                    sigmax = max(myvar{c}(timeframes,2));
                    posInds = find(myvar{c}(timeframes,2)>0);
                    sigmin = min(myvar{c}(posInds,2));
                end
                normvar(:,plotInd) = (myvar{c}(timeframes,2)-sigmin)/(sigmax-sigmin);
                normvar(:,plotInd+1) = (myvar{c+1}(timeframes,2)-sigmin)/(sigmax-sigmin);
            end
        else
            for c=plotind
                if get(handles.checkbox_commonNorm,'Value')
                    sigmax = str2num(get(handles.edit_maxI,'String'));
                    sigmin = str2num(get(handles.edit_minI,'String'));
                else
                    sigmax = max(myvar{c}(timeframes,2));
                    posInds = find(myvar{c}(timeframes,2)>0);
                    sigmin = min(myvar{c}(posInds,2));
                end
                normvar(:,plotInd) = (myvar{c}(timeframes,2)-sigmin)/(sigmax-sigmin);
                negInds = find(normvar(:,plotInd)<0);
                normvar(negInds,plotInd) = 0;
                plotInd=plotInd+1;
            end
        end
        
    case 3
        if pairdata==1
            for c=1:2:length(plotind)
                if get(handles.checkbox_commonNorm,'Value')
                    sigmin = str2num(get(handles.edit_minI,'String'));
                else
                    sigmax = max(myvar{c}(timeframes,2));
                end
                normvar(:,plotInd) = (myvar{c}(timeframes,2)-sigmin);
                normvar(:,plotInd) = (myvar{c+1}(timeframes,2)-sigmin);
            end
        else
            for c=plotind
                if get(handles.checkbox_commonNorm,'Value')
                    sigmin = str2num(get(handles.edit_minI,'String'));
                else
                    posInds = find(myvar{c}(timeframes,2)>0);
                    sigmin = min(myvar{c}(posInds,2));

                end
                normvar(:,plotInd) = (myvar{c}(timeframes,2)-sigmin);
                negInds = find(normvar(:,plotInd)<0);
                normvar(negInds,plotInd) = 0;
                plotInd=plotInd+1;
            end
        end
        
        
    case 4
        if pairdata==1
            for c=1:2:length(plotind)
                if get(handles.checkbox_commonNorm,'Value')
                    sigmin = str2num(get(handles.edit_minI,'String'));
                else
                    posInds = find(myvar{c}(timeframes,2)>0);
                    sigmin = min(myvar{c}(posInds,2));
                end
                normvar(:,plotInd) = (myvar{c}(timeframes,2)-sigmin);
                normvar(:,plotInd) = (myvar{c+1}(timeframes,2)-sigmin);
            end
        else
            for c=plotind
                
                if get(handles.checkbox_commonNorm,'Value')
                    sigmin = str2num(get(handles.edit_minI,'String'));
                else
                    posInds = find(myvar{c}(timeframes,2)>0);
                    sigmin = min(myvar{c}(posInds,2));
                end
                normvar(:,plotInd) = (myvar{c}(timeframes,2)-sigmin);
                negInds = find(normvar(:,plotInd)<0);
                normvar(negInds,plotInd) = 0;
                plotInd=plotInd+1;
            end
        end
end

kgrouping = kmeans(cellpath{end}(plotind,:),str2num(get(handles.edit_clusterno,'String')));

for c=1:size(normvar,2)
    if get(handles.popupmenu_plottype,'Value')==2;
        normvar(:,c) = smooth(normvar(:,c),str2num(get(handles.edit_framerange,'String')));
    end
    
    rankvalue(c,1) = c;
    [pks,locs] = findpeaks(normvar(:,c),'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',10,'THRESHOLD',0.7);
    rankvalue(c,2) = length(locs);
    rankvalue(c,3) = trapz(normvar(:,c));
    %rankvalue(c,4) = kgrouping(c);

end

neworder = sortrows(rankvalue,[3 2]);

load MyColormaps2
switch get(handles.popupmenu_colormap,'Value')
    case 1
        currentcmap = mycmap1;
    case 2
        currentcmap = mycmap2;
    case 3
        currentcmap = mycmap3;
    case 4
        currentcmap = mycmap4;
end

ind=1;
for t=timeframes
    ctime=round(myvar{plotind(1)}(t,1));
    if ~isempty(find(mod(ctime,60)==[0:6],1))
        plottime{ind,1} = num2str(round(ctime/60));
    else
        plottime{ind,1} = ' ';
    end
    ind=ind+1;
end

if ~get(handles.checkbox_rank,'Value')
    hmo = HeatMap((normvar)','Symmetric','false','ColumnLabels',plottime,...
        'RowLabels',plotind,'ColumnLabelsLocation','bottom','Colormap',currentcmap,...
        'ColumnLabelsRotate',0,'DisplayRange',str2num(get(handles.edit_maxnormI,'String')));
else
    %figure(),scatter(cellpath{end}(plotind,1),cellpath{end}(plotind,2),100,kgrouping,'filled');
    
    %for i=1:length(plotind)
    %    text(cellpath{end}(plotind(i),1)+15,cellpath{end}(plotind(i),2)+15,num2str(plotind(i)));hold on;
    %end
    hmo = HeatMap((normvar(:,neworder(:,1)))','Symmetric','false','ColumnLabels',plottime,...
        'RowLabels',plotind(neworder(:,1)),'ColumnLabelsLocation','bottom','Colormap',currentcmap,...
        'ColumnLabelsRotate',0,'DisplayRange',str2num(get(handles.edit_maxnormI,'String')));
    
end
addYLabel(hmo,'Cell No.');
addXLabel(hmo,'Time(hours)');

%cgo = clustergram((normvar'),'Symmetric','true','ColumnLabels',(handles.timestep),'RowLabels',plotind);
%get(cgo)
%datatypeInd = get(handles.popupmenu_datatype,'Value');
%availdatatypes = get(handles.popupmenu_datatype,'String');
%addTitle(hmo,['Cell population :' availdatatypes{datatypeInd}]);

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
