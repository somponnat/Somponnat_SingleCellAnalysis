function varargout = SignalClustering(varargin)
% SIGNALCLUSTERING MATLAB code for SignalClustering.fig
%      SIGNALCLUSTERING, by itself, creates a new SIGNALCLUSTERING or raises the existing
%      singleton*.
%
%      H = SIGNALCLUSTERING returns the handle to a new SIGNALCLUSTERING or the handle to
%      the existing singleton*.
%
%      SIGNALCLUSTERING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIGNALCLUSTERING.M with the given input arguments.
%
%      SIGNALCLUSTERING('Property','Value',...) creates a new SIGNALCLUSTERING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SignalClustering_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SignalClustering_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SignalClustering

% Last Modified by GUIDE v2.5 02-Aug-2013 09:19:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SignalClustering_OpeningFcn, ...
                   'gui_OutputFcn',  @SignalClustering_OutputFcn, ...
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


% --- Executes just before SignalClustering is made visible.
function SignalClustering_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SignalClustering (see VARARGIN)

% Choose default command line output for SignalClustering
handles.output = hObject;
handles.alldata = [];
handles.originData = [];
handles.cellFate = [];
handles.groupNo = [];
handles.timestamp = [];
handles.score = [];
handles.selectedcellIndices = [];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SignalClustering wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SignalClustering_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_posCtrl.
function popupmenu_posCtrl_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_posCtrl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_posCtrl contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_posCtrl
if isempty(handles.originData)
    set(handles.edit_commu,'String','No data. Please first initialize dataset.');
    return;
else
    axes(handles.axes_individual);
    plot(handles.timestamp,nanmean(handles.originData(handles.groupNo==get(hObject,'Value'),:),1),'g'); hold on;
    plot(handles.timestamp,nanmean(handles.originData(handles.groupNo==get(handles.popupmenu_negCtrl,'Value'),:),1),'r'); 
    plot(handles.timestamp,nanmean(handles.originData(handles.groupNo==get(handles.popupmenu_group,'Value'),:),1),'b');hold off;
    
    table_data = get(handles.uitable_params,'Data');
    for i=[1:3 5:size(table_data,1)]
        table_data{i,3} = nanmean(handles.alldata(handles.groupNo==get(hObject,'Value'),i));
    end
    table_data{4,3} = [];
    set(handles.uitable_params,'Data',table_data);
end

if ~isempty(handles.score)
    axes(handles.axes_posctrl);
    countInd=1;
    binSize=[];
    for j=handles.newList
        binSize(countInd) = numel(handles.binning{get(handles.popupmenu_posCtrl,'Value'),j});
        countInd=countInd+1;
    end
    bar(1:str2num(get(handles.edit_clusterno,'String')),binSize);
    set(handles.axes_posctrl,'XTickLabel',[]);
end

% --- Executes during object creation, after setting all properties.
function popupmenu_posCtrl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_posCtrl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_y_pca.
function popupmenu_y_pca_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_y_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_y_pca contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_y_pca

[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;
guidata(hObject, handles); 
% --- Executes during object creation, after setting all properties.
function popupmenu_y_pca_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_y_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
[filename,PathName,FilterIndex] = uigetfile('*.nd', 'Choose metamorph ND file','C:\computation\02-03-2013\02032013-r1.nd');
if FilterIndex~=0
    set(handles.edit_ndfilename,'String',filename);
    handles.ndfilename = filename;
    handles.ndpathname = PathName;
    set(handles.edit_sourceF,'String',PathName); 
end
guidata(hObject, handles);  

% --- Executes on button press in pushbutton_initialize.
function pushbutton_initialize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_initialize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.ndpathname)
    set(handles.edit_commu,'String','Please first choose ND file');
    return;
end

searchInd(1,:)= str2num(get(handles.edit_selectedRows,'String'));
searchInd(2,:)= str2num(get(handles.edit_selectedCols,'String'));
searchInd(3,:)= str2num(get(handles.edit_selectedFields,'String'));

phenotypeList = {'Dead(-4)';'Dead,Divided 3 times(-3)';'Dead,Divided twice(-2)';'Dead,Divided once(-1)';'Quiescent(0)';'Divided once(1)';'Divided twice(2)';'Divided 3 times(3)'};

outputsignalNo = str2num(get(handles.edit_outputsignalno,'String'));
alldata=[];
groupNo=[];
originData=[];
cellFate=[];

for i = 1:size(searchInd,2)
        switch mod(i,3)
            case 1
                set(handles.edit_commu,'String','Processing.');
            case 2
                set(handles.edit_commu,'String','Processing..');
            otherwise
                set(handles.edit_commu,'String','Processing...');
        end
        pause(0.01);
        
        row = searchInd(1,i);
        col = searchInd(2,i);
        field = searchInd(3,i);
        H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
        param_name = ['/field' num2str(field)  '/clusterparams' num2str(outputsignalNo)];
        
        signal_name = ['/field' num2str(field)  '/outputsignal' num2str(outputsignalNo)];
        timestamp_name = ['/field' num2str(field) '/timestamp' num2str(outputsignalNo)];

        fid = H5F.open(fullfile(handles.ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
        if H5L.exists(fid,param_name,'H5P_DEFAULT')
            H5F.close(fid);
            paraminfo = h5info(fullfile(handles.ndpathname,H5filename), param_name);
            startind = double([1 1]);
            countind = [paraminfo.Dataspace.Size(1) paraminfo.Dataspace.Size(2)];
            param_mat = double(h5read(fullfile(handles.ndpathname,H5filename),param_name,startind, countind));
            alldata = [alldata;param_mat];
            groupNo = [groupNo;i*ones(size(param_mat,1),1)];
            
            signalinfo = h5info(fullfile(handles.ndpathname,H5filename), signal_name);
            %sisterListinfo = h5info(fullfile(handles.SourceF,H5filename), sisterList_name);
            %sisterList = h5read(fullfile(handles.SourceF,H5filename),sisterList_name,[1 1 1],...
            %                 [sisterListinfo.Dataspace.Size(1) sisterListinfo.Dataspace.Size(2) sisterListinfo.Dataspace.Size(3)]);
            startind = double([1 1 2]);
            countind = [signalinfo.Dataspace.Size(1) signalinfo.Dataspace.Size(2) 1];
            signal = permute(h5read(fullfile(handles.ndpathname,H5filename),signal_name,startind, countind),[2 1 3]);

            for c_cell=1:size(param_mat,1)
                originData = [originData;signal(:,param_mat(c_cell,4))'];
                cellFate =  [cellFate;param_mat(c_cell,5)];
            end
            if i==1
                timestamp = h5read(fullfile(handles.ndpathname,H5filename),timestamp_name);
                table_data = cell(size(param_mat,2),5);
                for j = 1:size(param_mat,2)
                    table_data{j,1} = h5readatt(fullfile(handles.ndpathname,H5filename),param_name,['param' num2str(j)]);
                end
            end
            clear param_mat;
            
            
        else
            H5F.close(fid);
        end
end

[R C] = find(originData==0);
for i=1:length(R)
    originData(R(i),C(i)) = NaN;
end

set(handles.togglebutton_showselectedpolygon,'Value',0);
handles.plottype = 1;
set(handles.togglebutton_scatter2D,'Value',1);
handles.selectedcellIndices = [];
handles.selectedPolygon = [];
handles.alldata = alldata;
handles.originData = originData;
handles.phenotypeList = phenotypeList;

Phenolabel{1} = 'All';
PhenotypeID = unique(cellFate);
for i=1:length(PhenotypeID)
    Phenolabel{i+1} = [phenotypeList{PhenotypeID(i)+5}];
end

Groupnolabel{1} = 'All';
groupNoID = unique(groupNo);
for i=1:length(groupNoID)
    Groupnolabel{i+1} = ['r' num2str(searchInd(1,groupNoID(i))) 'c' num2str(searchInd(2,groupNoID(i)))];
end

Clusterlabel{1} = 'All';
ClusterNo = str2num(get(handles.edit_clusterno,'String'));
for i=1:ClusterNo
    Clusterlabel{i+1} = [num2str(i)];
end

set(handles.popupmenu_phenotype,'String',Phenolabel);
set(handles.popupmenu_wellposition,'String',Groupnolabel);
set(handles.popupmenu_cluster,'String',Clusterlabel);

handles.cellFate = cellFate;
choiceInput = sort(unique(cellFate))';
set(handles.edit_phetype_list,'String',num2str(choiceInput));

handles.groupNo = groupNo;
choiceInput = sort(unique(groupNo))';
set(handles.edit_position_list,'String',num2str(choiceInput));

choiceInput = 1:str2num(get(handles.edit_clusterno,'String'));
set(handles.edit_cluster_list,'String',num2str(choiceInput));

handles.timestamp = timestamp;
handles.plotInd = 1:size(alldata,1);
handles.grayInd = [];
guidata(hObject, handles);
handles = guidata(hObject);
axes(handles.axes_individual);
plot(handles.timestamp,nanmean(handles.originData(handles.groupNo==get(handles.popupmenu_posCtrl,'Value'),:),1),'g'); hold on;
plot(handles.timestamp,nanmean(handles.originData(handles.groupNo==get(handles.popupmenu_negCtrl,'Value'),:),1),'r');
plot(handles.timestamp,nanmean(handles.originData(handles.groupNo==get(handles.popupmenu_group,'Value'),:),1),'b');hold off;

for i=[1:3 5:size(table_data,1)]
    table_data{i,3} = nanmean(handles.alldata(handles.groupNo==get(handles.popupmenu_posCtrl,'Value'),i));
end
table_data{4,3} = [];
for i=[1:3 5:size(table_data,1)]
    table_data{i,4} = nanmean(handles.alldata(handles.groupNo==get(handles.popupmenu_negCtrl,'Value'),i));
end
table_data{4,4} = [];
for i=[1:3 5:size(table_data,1)]
    table_data{i,5} = nanmean(handles.alldata(handles.groupNo==get(handles.popupmenu_group,'Value'),i));
end
table_data{4,5} = [];

set(handles.uitable_params,'Data',table_data);
set(handles.edit_commu,'String','Finished initializing parameters');
guidata(hObject, handles);

selectedparams = str2num(get(handles.edit_selectedparams,'String'));
table_data = get(handles.uitable_params,'Data');
param_names = [];

for i=1:size(table_data,1)
    param_names{i} = table_data{i,1};
end
Ind=1;
for i=selectedparams
    selectedparams_name{Ind} = param_names{i};
    Ind = Ind+1;
end
set(handles.popupmenu_selectedparams,'String',selectedparams_name);




function edit_sourceF_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sourceF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sourceF as text
%        str2double(get(hObject,'String')) returns contents of edit_sourceF as a double
handles.ndpathname = get(hObject,'String');
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



function edit_selectedRows_Callback(hObject, eventdata, handles)
% hObject    handle to edit_selectedRows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_selectedRows as text
%        str2double(get(hObject,'String')) returns contents of edit_selectedRows as a double


% --- Executes during object creation, after setting all properties.
function edit_selectedRows_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_selectedRows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_selectedCols_Callback(hObject, eventdata, handles)
% hObject    handle to edit_selectedCols (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_selectedCols as text
%        str2double(get(hObject,'String')) returns contents of edit_selectedCols as a double


% --- Executes during object creation, after setting all properties.
function edit_selectedCols_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_selectedCols (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_group.
function popupmenu_group_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_group contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_group
if isempty(handles.originData)
    set(handles.edit_commu,'String','No data. Please first initialize dataset.');
    return;
else
    axes(handles.axes_individual);
    plot(handles.timestamp,nanmean(handles.originData(handles.groupNo==get(handles.popupmenu_negCtrl,'Value'),:),1),'r'); hold on;
    plot(handles.timestamp,nanmean(handles.originData(handles.groupNo==get(handles.popupmenu_posCtrl,'Value'),:),1),'g');
    plot(handles.timestamp,nanmean(handles.originData(handles.groupNo==get(hObject,'Value'),:),1),'b'); hold off;
    
    table_data = get(handles.uitable_params,'Data');
    for i=[1:3 5:size(table_data,1)]
        table_data{i,5} = nanmean(handles.alldata(handles.groupNo==get(hObject,'Value'),i));
    end
    table_data{4,5} = [];
    set(handles.uitable_params,'Data',table_data);
end

if ~isempty(handles.score)
    axes(handles.axes_group);
    countInd=1;
    binSize=[];
    for j=handles.newList
        binSize(countInd) = numel(handles.binning{get(handles.popupmenu_group,'Value'),j});
        countInd=countInd+1;
    end
    bar(1:str2num(get(handles.edit_clusterno,'String')),binSize);
    set(handles.axes_group,'XTickLabel',[]);
end

% --- Executes during object creation, after setting all properties.
function popupmenu_group_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_negCtrl.
function popupmenu_negCtrl_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_negCtrl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_negCtrl contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_negCtrl
if isempty(handles.originData)
    set(handles.edit_commu,'String','No data. Please first initialize dataset.');
    return;
else
    axes(handles.axes_individual);
    plot(handles.timestamp,nanmean(handles.originData(handles.groupNo==get(hObject,'Value'),:),1),'r'); hold on;
    plot(handles.timestamp,nanmean(handles.originData(handles.groupNo==get(handles.popupmenu_posCtrl,'Value'),:),1),'g');
    plot(handles.timestamp,nanmean(handles.originData(handles.groupNo==get(handles.popupmenu_group,'Value'),:),1),'b');hold off;
    
    table_data = get(handles.uitable_params,'Data');
    for i=[1:3 5:size(table_data,1)]
        table_data{i,4} = nanmean(handles.alldata(handles.groupNo==get(hObject,'Value'),i));
    end
    table_data{4,4} = [];
    set(handles.uitable_params,'Data',table_data);
end

if ~isempty(handles.score)
    axes(handles.axes_negctrl);
    countInd=1;
    binSize=[];
    for j=handles.newList
        binSize(countInd) = numel(handles.binning{get(handles.popupmenu_negCtrl,'Value'),j});
        countInd=countInd+1;
    end
    bar(1:str2num(get(handles.edit_clusterno,'String')),binSize);
    set(handles.axes_negctrl,'XTickLabel',[]);
end

% --- Executes during object creation, after setting all properties.
function popupmenu_negCtrl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_negCtrl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function [notp stagePos stageName waveName] = readndfile(pathname,filename)
% Search for number of string matches per line.  
notp=-1;
stagePos = [];
stageName = [];
waveName = [];
currentF = pwd;


if exist(fullfile(pathname,filename),'file')
    fid = fopen(fullfile(pathname,filename));
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
            stagename  = regexp(tline, '(?<="Stage\d+", ").+(?=")', 'match');
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


% --- Executes on button press in pushbutton_loadbyndfile.
function pushbutton_loadbyndfile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadbyndfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

stageName = [];
searchInd(1,:)= str2num(get(handles.edit_selectedRows,'String'));
searchInd(2,:)= str2num(get(handles.edit_selectedCols,'String'));
searchInd(3,:)= str2num(get(handles.edit_selectedFields,'String'));
for i=1:size(searchInd,2)
    stageName{i} = ['r' num2str(searchInd(1,i)) 'c' num2str(searchInd(2,i)) 'f' num2str(searchInd(3,i))];
end


set(handles.popupmenu_posCtrl,'String',stageName);
set(handles.popupmenu_posCtrl,'Value',1);
set(handles.popupmenu_negCtrl,'String',stageName);
set(handles.popupmenu_negCtrl,'Value',1);
set(handles.popupmenu_group,'String',stageName);
set(handles.popupmenu_group,'Value',1);


guidata(hObject, handles);  
set(handles.edit_commu,'String',['Assigned source folder to ' handles.ndpathname]);
% --- Executes on button press in pushbutton_selectPoints.
function pushbutton_selectPoints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit_commu,'String','Draw polygon to choose cells. Click last point with Right mouse button.');
searchInd(1,:)= str2num(get(handles.edit_selectedRows,'String'));
searchInd(2,:)= str2num(get(handles.edit_selectedCols,'String'));
searchInd(3,:)= str2num(get(handles.edit_selectedFields,'String'));
old_xlim = get(handles.axes_pca,'XLim');
old_ylim = get(handles.axes_pca,'YLim');
warning off;
[plot_h,plotInd,grayInd] = plotPCA(0,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;

set(handles.axes_pca,'XLim',old_xlim);
set(handles.axes_pca,'YLim',old_ylim);
guidata(hObject, handles);  
handles = guidata(hObject);  

hold on
% Initially, the list of points is empty.
xy = [];
n = 0;
% Loop, picking up the points.

but = 1;
while but == 1
    [xi,yi,but] = ginput(1);
    n = n+1;
    xy(n,:) = [xi yi];
    if n==1
        plot(xi,yi,'k-');
    else
        plot(xy([n-1 n],1),xy([n-1 n],2),'k-');
    end
end

plot(xy([1 n],1),xy([1 n],2),'k-');

for i=1:length(handles.plot_h)
    c_X = get(handles.plot_h,'XData')';
    c_Y = get(handles.plot_h,'YData')';
    IN= inpolygon(c_X,c_Y,xy(:,1),xy(:,2));
end
hold off


insideInd = find(IN==1);
mycolor = jet(size(searchInd,2));
axes(handles.axes_individual);
for i=1:length(insideInd)
    plot(handles.timestamp(handles.originData(handles.plotInd(insideInd(i)),:)~=0),handles.originData(handles.plotInd(insideInd(i)),handles.originData(handles.plotInd(insideInd(i)),:)~=0),'Color',mycolor(handles.groupNo(handles.plotInd(insideInd(i))),:));hold on;
end
hold off;
drawnow;
warning on;

handles.selectedcellIndices = handles.plotInd(insideInd);
handles.selectedPolygon = xy;
guidata(hObject, handles);  
handles = guidata(hObject);

table_data = get(handles.uitable_params,'Data');
for i=1:size(table_data,1)
    table_data{i,2} = nanmean(handles.alldata(handles.plotInd(insideInd),i));
end
set(handles.uitable_params,'Data',table_data);
[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;

hold on;
plot(xy([1:n 1],1),xy([1:n 1],2),'k-');hold off;
set(handles.axes_pca,'XLim',old_xlim );
set(handles.axes_pca,'YLim',old_ylim );

noCluster = str2num(get(handles.edit_clusterno,'String'));
binSize=[];
for j=1:noCluster
    
    binSize(j) = numel(find(handles.T(handles.plotInd(insideInd))==j));
end

axes(handles.axes_selected);
bar(1:str2num(get(handles.edit_clusterno,'String')),binSize);
set(handles.axes_selected,'XTickLabel',[]);
set(handles.togglebutton_showselectedpolygon,'Value',1);
set(handles.edit_commu,'String',['Cells in polygon: ' num2str(length(handles.plotInd(insideInd))) ' of ' num2str(length(handles.plotInd))]);
guidata(hObject, handles);  

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_removepoints.
function pushbutton_removepoints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removepoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [plotInd,grayInd] = determinePlotInd(handles)
allInd = 1:size(handles.alldata,1);
PhenoInd = str2num(get(handles.edit_phetype_list,'String'));
WellInd = str2num(get(handles.edit_position_list,'String'));
CInd = str2num(get(handles.edit_cluster_list,'String'));
SInd = get(handles.popupmenu_selection,'Value');
if SInd == 2
    selectedCells = handles.selectedcellIndices;
else
    selectedCells = allInd;
end
selectedPheno = [];
for i=1:length(PhenoInd)
    selectedPheno = [selectedPheno;find(handles.cellFate == PhenoInd(i))];
end
selectedPheno = sort(selectedPheno);

selectedWell = [];
for i=1:length(WellInd)
    selectedWell = [selectedWell;find(handles.groupNo == WellInd(i))];
end
selectedWell = sort(selectedWell);

selectedC = [];
for i=1:length(CInd)
    selectedC = [selectedC;find(handles.T == CInd(i))];
end
selectedC = sort(selectedC);

plotInd = sort(intersect(intersect(intersect(selectedPheno,selectedWell),selectedC),selectedCells));
grayInd = sort(setdiff(1:size(handles.alldata,1),plotInd));
set(handles.edit_totalcellcount,'String',num2str(length(plotInd)));
set(handles.edit_commu,'String',['Total cell selected: ' num2str(length(plotInd))]);

function [plot_h,plotInd,grayInd]= plotPCA(plottype,handles,selectedcellInd)
[plotInd,grayInd] = determinePlotInd(handles);
x = get(handles.popupmenu_x_pca,'Value');
y = get(handles.popupmenu_y_pca,'Value');
z = get(handles.popupmenu_z_pca,'Value');
c = get(handles.popupmenu_c_pca,'Value');
p = get(handles.popupmenu_selectedparams,'Value');
selectedparams = str2num(get(handles.edit_selectedparams,'String'));

if get(handles.togglebutton_orthorotate,'Value')==0
    myplotdata = handles.score;
else
    myplotdata = handles.rotated_score;
end



msize = str2num(get(handles.edit_markersize,'String'));
mtype = get(handles.edit_markertype,'String');

plot_h=[];
searchInd(1,:)= str2num(get(handles.edit_selectedRows,'String'));
searchInd(2,:)= str2num(get(handles.edit_selectedCols,'String'));
searchInd(3,:)= str2num(get(handles.edit_selectedFields,'String'));

if ~isempty(myplotdata)
    axes(handles.axes_pca);
    
    switch plottype
        case 1 % scatter-2D
            if ~isempty(selectedcellInd) && get(handles.togglebutton_showselectedpolygon,'Value')==1;
                scatter(myplotdata(intersect(plotInd,selectedcellInd),x),myplotdata(intersect(plotInd,selectedcellInd),y),msize,'k','o','fill');hold on;
            end
            if ~isempty(grayInd)
                scatter(myplotdata(grayInd,x),myplotdata(grayInd,y),msize,[0.7 0.7 0.7],'x');hold on;
            end
            if ~isempty(plotInd)
                switch c
                    case 1 % well position
                        c_groupNo = unique(handles.groupNo(plotInd));
                        for i=1:length(c_groupNo)
                            myLegend{i} = ['r' num2str(searchInd(1,c_groupNo(i))) 'c' num2str(searchInd(2,c_groupNo(i)))];
                        end
                        
                        plot_h = gscatter(myplotdata(plotInd,x),myplotdata(plotInd,y),nominal(handles.groupNo(plotInd),myLegend),jet(length(unique(handles.groupNo(plotInd)))),mtype,msize-6,'on');
                        set(handles.togglebutton_legendLogic,'Value',1);
                        
                    case 2 % phenotype
                        c_phenotype = unique(handles.cellFate(plotInd));
                        mymarkertype = [];
                        mymarkercolor = [];
                        for i=1:length(c_phenotype)
                            myLegend{i} = handles.phenotypeList{c_phenotype(i)+5};
                            switch c_phenotype(i)
                                case -4
                                    mymarkertype = [mymarkertype 'v'];
                                case {-3,-2,-1}
                                    mymarkertype = [mymarkertype 'v'];
                                case 0
                                    mymarkertype = [mymarkertype 'o'];
                                case {1,2,3}
                                    mymarkertype = [mymarkertype 'x'];
                            end
                            
                            switch c_phenotype(i)
                                case -4
                                    mymarkercolor = [mymarkercolor 'k'];
                                case {-3,3}
                                    mymarkercolor = [mymarkercolor 'b'];
                                case {-2,2}
                                    mymarkercolor = [mymarkercolor 'g'];
                                case {-1,1}
                                    mymarkercolor = [mymarkercolor 'r'];
                                case 0
                                    mymarkercolor = [mymarkercolor 'c'];
                            end
                        end
                        plot_h = gscatter(myplotdata(plotInd,x),myplotdata(plotInd,y),nominal(handles.cellFate(plotInd),myLegend),mymarkercolor,mymarkertype,msize-6,'on');
                        set(handles.togglebutton_legendLogic,'Value',1);
                        %legend(myLegend);
                    case 3 % cluster
                        plot_h = gscatter(myplotdata(plotInd,x),myplotdata(plotInd,y),handles.T(plotInd),hsv(length(unique(handles.T(plotInd)))),mtype,msize-6,'on');
                        set(handles.togglebutton_legendLogic,'Value',1);
                    case 4 % plot params
                        chosen_param = selectedparams(p);
                        alldataset = handles.alldata(:,chosen_param);
                        myx = [min(alldataset):(max(alldataset)-min(alldataset))/31:max(alldataset)];
                        [~,bin] = histc(handles.alldata(plotInd,chosen_param),myx);
                        mycolor = hot(32);
                        for i=1:length(bin)
                            cellcolor(i,:) = mycolor(bin(i),:);
                        end
                        plot_h = scatter(myplotdata(plotInd,x),myplotdata(plotInd,y),msize,cellcolor,mtype);
                        set(handles.togglebutton_legendLogic,'Value',0);
                end
                hold off;
                colorbar('off');
            end
        case 2
            if ~isempty(selectedcellInd) && get(handles.togglebutton_showselectedpolygon,'Value')==1;
                scatter3(myplotdata(intersect(plotInd,selectedcellInd),x),myplotdata(intersect(plotInd,selectedcellInd),y),myplotdata(intersect(plotInd,selectedcellInd),z),msize,'k','o','fill');hold on;
            end
            if ~isempty(grayInd)
                scatter3(myplotdata(grayInd,x),myplotdata(grayInd,y),myplotdata(grayInd,z),msize,[0.7 0.7 0.7],'x');hold on;
            end
            if ~isempty(plotInd)
                switch c
                    case 1 % well position
                        set(handles.edit_commu,'String','Plot colors show well positions.');
                        cellcolor = [];
                        mycolor = jet(size(searchInd,2));
                        mydata = handles.groupNo(plotInd);
                        for i=1:length(mydata)
                            cellcolor(i,:) = mycolor(mydata(i),:);
                        end
                        plot_h=scatter3(myplotdata(plotInd,x),myplotdata(plotInd,y),myplotdata(plotInd,z),msize,cellcolor,mtype);
                        set(handles.togglebutton_legendLogic,'Value',0);
                    case 2 % phenotype
                        set(handles.edit_commu,'String','Plot colors show cell decision.');
                        
                        plot_h=scatter3(myplotdata(intersect(plotInd,find(handles.cellFate==-4)),x),myplotdata(intersect(plotInd,find(handles.cellFate==-4)),y),myplotdata(intersect(plotInd,find(handles.cellFate==-4)),z),msize,'v','k');hold on;
                        plot_h=scatter3(myplotdata(intersect(plotInd,find(handles.cellFate==-3)),x),myplotdata(intersect(plotInd,find(handles.cellFate==-3)),y),myplotdata(intersect(plotInd,find(handles.cellFate==-3)),z),msize,'v','b');
                        plot_h=scatter3(myplotdata(intersect(plotInd,find(handles.cellFate==-2)),x),myplotdata(intersect(plotInd,find(handles.cellFate==-2)),y),myplotdata(intersect(plotInd,find(handles.cellFate==-2)),z),msize,'v','g');
                        plot_h=scatter3(myplotdata(intersect(plotInd,find(handles.cellFate==-1)),x),myplotdata(intersect(plotInd,find(handles.cellFate==-1)),y),myplotdata(intersect(plotInd,find(handles.cellFate==-1)),z),msize,'v','r');
                        plot_h=scatter3(myplotdata(intersect(plotInd,find(handles.cellFate==0)),x),myplotdata(intersect(plotInd,find(handles.cellFate==0)),y),myplotdata(intersect(plotInd,find(handles.cellFate==0)),z),msize,'o','c');
                        plot_h=scatter3(myplotdata(intersect(plotInd,find(handles.cellFate==1)),x),myplotdata(intersect(plotInd,find(handles.cellFate==1)),y),myplotdata(intersect(plotInd,find(handles.cellFate==1)),z),msize,'x','r');
                        plot_h=scatter3(myplotdata(intersect(plotInd,find(handles.cellFate==2)),x),myplotdata(intersect(plotInd,find(handles.cellFate==2)),y),myplotdata(intersect(plotInd,find(handles.cellFate==2)),z),msize,'x','g');
                        plot_h=scatter3(myplotdata(intersect(plotInd,find(handles.cellFate==3)),x),myplotdata(intersect(plotInd,find(handles.cellFate==3)),y),myplotdata(intersect(plotInd,find(handles.cellFate==3)),z),msize,'x','b');hold off;
                        
                        set(handles.togglebutton_legendLogic,'Value',0);
                    case 3 % cluster
                        set(handles.edit_commu,'String','Plot colors show cluster number.');
                        cellcolor = [];
                        mycolor = hsv(length(unique(handles.T(plotInd))));
                        mydata = handles.T(plotInd);
                        for i=1:length(mydata)
                            cellcolor(i,:) = mycolor(mydata(i),:);
                        end
                        plot_h=scatter3(myplotdata(plotInd,x),myplotdata(plotInd,y),myplotdata(plotInd,z),msize,cellcolor,mtype);
                        set(handles.togglebutton_legendLogic,'Value',0);
                    case 4 % plot params
                        set(handles.edit_commu,'String','Plot colors show intensity of the selected parameter.');
                        
                        chosen_param = selectedparams(p);
                        alldataset = handles.alldata(:,chosen_param);
                        myx = [min(alldataset):(max(alldataset)-min(alldataset))/31:max(alldataset)];
                        [~,bin] = histc(handles.alldata(plotInd,chosen_param),myx);
                        mycolor = hot(32);
                        for i=1:length(bin)
                            cellcolor(i,:) = mycolor(bin(i),:);
                        end
                        plot_h=scatter3(myplotdata(plotInd,x),myplotdata(plotInd,y),myplotdata(plotInd,z),msize,cellcolor,mtype);
                        set(handles.togglebutton_legendLogic,'Value',0);
                end
                hold off;
                xnames = get(handles.popupmenu_x_pca,'String');
                ynames = get(handles.popupmenu_y_pca,'String');
                znames = get(handles.popupmenu_z_pca,'String');
                xlabel(xnames{get(handles.popupmenu_x_pca,'Value')});
                ylabel(ynames{get(handles.popupmenu_y_pca,'Value')});
                zlabel(znames{get(handles.popupmenu_z_pca,'Value')});
            end
        case 3
            if ~isempty(plotInd)
                set(handles.togglebutton_legendLogic,'Value',0);
                x=myplotdata(plotInd,x);
                y=myplotdata(plotInd,y);
                out = scatplot(x,y,'circles',str2num(get(handles.edit_contourcirclesize,'String')),str2num(get(handles.edit_meshsize,'String')),8,2,msize);
                view(2);
                xlabel([]);
                ylabel([]);
                zlabel([]);
            end

        case 0 
            set(handles.togglebutton_legendLogic,'Value',0);
            switch c
                case 1 % well position
                    
                    cellcolor = [];
                    mycolor = jet(size(searchInd,2));
                    mydata = handles.groupNo(plotInd);
                    for i=1:length(mydata)
                        cellcolor(i,:) = mycolor(mydata(i),:);
                    end
                    plot_h = scatter(myplotdata(plotInd,x),myplotdata(plotInd,y),msize,cellcolor,mtype);
                                        
                case 2 % phenotype

                    cellcolor = [];
                    mycolor = lines(7);
                    mydata = handles.cellFate(plotInd);
                    for i=1:length(mydata)
                        cellcolor(i,:) = mycolor(mydata(i)+5,:);
                    end
                    plot_h = scatter(myplotdata(plotInd,x),myplotdata(plotInd,y),msize,cellcolor,mtype);
                    
                case 3 % cluster
                    cellcolor = [];
                    mycolor = hsv(str2num(get(handles.edit_clusterno,'String')));
                    mydata = handles.T(plotInd);
                    for i=1:length(mydata)
                        cellcolor(i,:) = mycolor(mydata(i),:);
                    end
                    plot_h = scatter(myplotdata(plotInd,x),myplotdata(plotInd,y),msize,cellcolor,mtype);
                case 4 % plot params
                    chosen_param = selectedparams(p);
                    alldataset = handles.alldata(:,chosen_param);
                    myx = [min(alldataset):(max(alldataset)-min(alldataset))/31:max(alldataset)];
                    [~,bin] = histc(handles.alldata(plotInd,chosen_param),myx);
                    mycolor = hot(32);
                    for i=1:length(bin)
                        cellcolor(i,:) = mycolor(bin(i),:);
                    end
                    plot_h = scatter(myplotdata(plotInd,x),myplotdata(plotInd,y),msize,cellcolor,mtype);
            end
    end
    hold off;
else
    set(handles.edit_commu,'String','Must first process data');
end


% --- Executes on selection change in popupmenu_x_pca.
function popupmenu_x_pca_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_x_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_x_pca contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_x_pca

[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;
guidata(hObject, handles); 

% --- Executes during object creation, after setting all properties.
function popupmenu_x_pca_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_x_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_group2.
function popupmenu_group2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_group2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_group2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_group2


% --- Executes during object creation, after setting all properties.
function popupmenu_group2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_group2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_group3.
function popupmenu_group3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_group3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_group3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_group3


% --- Executes during object creation, after setting all properties.
function popupmenu_group3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_group3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_selectedparams.
function popupmenu_selectedparams_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_selectedparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_selectedparams contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_selectedparams
set(handles.popupmenu_c_pca,'Value',4);
[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function popupmenu_selectedparams_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_selectedparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu16.
function popupmenu16_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu16 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu16


% --- Executes during object creation, after setting all properties.
function popupmenu16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_c_pca.
function popupmenu_c_pca_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_c_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_c_pca contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_c_pca

[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;

guidata(hObject, handles); 

% --- Executes during object creation, after setting all properties.
function popupmenu_c_pca_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_c_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_addparams.
function pushbutton_addparams_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_removeparams.
function pushbutton_removeparams_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removeparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_showpcaloading.
function pushbutton_showpcaloading_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_showpcaloading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.togglebutton_orthorotate,'Value')==0
    myplotdata = handles.score;
    mycoefforth = handles.coefforth;
else
    myplotdata = handles.rotated_score;
    mycoefforth = handles.rotated_coefforth;
end



figure(1);
paramNo = str2num(get(handles.edit_selectedparams,'String'));
subplot(2,3,1);
biplot(mycoefforth(:,1:3),'scores',myplotdata(:,1:3),'varlabels',get(handles.popupmenu_selectedparams,'String'));
subplot(2,3,2:3);
pareto(handles.explained);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
s(1) = subplot(2,4,5);
barh(mycoefforth(:,1));title('PC1');
set(gca,'YLim',[0 length(paramNo)+1],'YTick',1:length(paramNo),'YTickLabel',get(handles.popupmenu_selectedparams,'String'));
s(2) = subplot(2,4,6);
barh(mycoefforth(:,2));title('PC2');
set(gca,'YTickLabel',[]);
s(3) = subplot(2,4,7);
barh(mycoefforth(:,3));title('PC3');
set(gca,'YTickLabel',[]);
s(4) = subplot(2,4,8);
barh(mycoefforth(:,4));title('PC4');
set(gca,'YTickLabel',[]);
linkaxes(s,'y');
% --- Executes on button press in pushbutton_cluster.
function pushbutton_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedparams = str2num(get(handles.edit_selectedparams,'String'));
noCluster = str2num(get(handles.edit_clusterno,'String'));
searchInd(1,:)= str2num(get(handles.edit_selectedRows,'String'));
searchInd(2,:)= str2num(get(handles.edit_selectedCols,'String'));
searchInd(3,:)= str2num(get(handles.edit_selectedFields,'String'));

w = 1./nanvar(handles.alldata(:,selectedparams));
%w(w==inf) = 1e-64;
[wcoeff,score,latent,tsquared,explained] = pca(handles.alldata(:,selectedparams),'VariableWeights',w);
[wcoeff_rotated,rotationMat] = rotatefactors(wcoeff);
coefforth = diag(sqrt(w))*wcoeff;
rotated_coefforth = diag(sqrt(w))*wcoeff_rotated;
rotated_score = zscore(handles.alldata(:,selectedparams))*rotated_coefforth;
set(handles.edit_commu,'String','Finished calculating PCA');
handles.explained = explained;
handles.score = score;
handles.rotated_score = rotated_score;
handles.coefforth = coefforth;
handles.rotated_coefforth = rotated_coefforth;

if get(handles.togglebutton_orthorotate,'Value')==0
    Z = linkage(score(:,1:3),'ward','euclidean');
else
    Z = linkage(rotated_score(:,1:3),'ward','euclidean');
end

        
T = cluster(Z,'maxclust',noCluster);        
binning = [];
for i=1:size(searchInd,2)
    binSize=[];
    for j=1:noCluster
        binning{i,j} = find(handles.groupNo==i & T==j);
        binSize(j) = numel(binning{i,j});
    end
end
clear refGroupList myList1 taken_myList1 n x temp sizeTaken newList

refGroup = get(handles.popupmenu_negCtrl,'Value');
refGroupList = T(handles.groupNo==refGroup);
uniqueCluster = unique(refGroupList);
temp = zeros(length(uniqueCluster),2);
for i=1:length(uniqueCluster)
    temp(i,1) = uniqueCluster(i);
    temp(i,2) = length(find(refGroupList==uniqueCluster(i)));
end

sizeTaken =length(uniqueCluster);
myList1 = flipud(sortrows(temp,2));
if size(myList1,1)<sizeTaken
    taken_myList1 = myList1(1:size(myList1,1),1)';
else
    taken_myList1 = myList1(1:sizeTaken,1)';
end

clear refGroupList taken_myList2  temp sizeTaken2 

refGroup = get(handles.popupmenu_posCtrl,'Value');

refGroupList = T(handles.groupNo==refGroup);
uniqueCluster = unique(refGroupList);
temp2 = zeros(length(uniqueCluster),2);
for i=1:length(uniqueCluster)
    temp2(i,1) = uniqueCluster(i);
    temp2(i,2) = length(find(refGroupList==uniqueCluster(i)));
end

sizeTaken2 =length(uniqueCluster);
myList2 = flipud(sortrows(temp2,2));


if size(myList1,1)<sizeTaken2
    taken_myList2 = myList2(1:size(myList2,1),1)';
else
    taken_myList2 = myList2(1:sizeTaken2,1)';
end
clear newList;
newList = zeros(1,noCluster);
FrontI = 1;
BackI = noCluster;
for i=1:max([length(taken_myList1) length(taken_myList2)])
    if i<=length(taken_myList1) && isempty(find(newList==taken_myList1(i),1))
        newList(FrontI) = taken_myList1(i);
        taken_myList2 = taken_myList2(taken_myList2~=taken_myList1(i));
        FrontI=FrontI+1;
    end
    
    if i<=length(taken_myList2) && isempty(find(newList==taken_myList2(i),1))
        newList(BackI) = taken_myList2(i);
        taken_myList1 = taken_myList1(taken_myList1~=taken_myList2(i));
        BackI=BackI-1;
    end
    
end
emptyInd = find(newList==0);
emptyValue = setdiff(1:noCluster,newList);
for i=1:length(emptyInd)
    newList(emptyInd(i)) = emptyValue(i);
end

axes(handles.axes_posctrl);
countInd=1;
binSize=[];
for j=newList
    binSize(countInd) = numel(binning{get(handles.popupmenu_posCtrl,'Value'),j});
    countInd=countInd+1;
end
bar(1:noCluster,binSize);
set(handles.axes_posctrl,'XTickLabel',[]);

axes(handles.axes_negctrl);
countInd=1;
binSize=[];
for j=newList
    binSize(countInd) = numel(binning{get(handles.popupmenu_negCtrl,'Value'),j});
    countInd=countInd+1;
end
bar(1:noCluster,binSize);
set(handles.axes_negctrl,'XTickLabel',[]);

axes(handles.axes_group);
countInd=1;
binSize=[];
for j=newList
    binSize(countInd) = numel(binning{get(handles.popupmenu_group,'Value'),j});
    countInd=countInd+1;
end
bar(1:noCluster,binSize);
set(handles.axes_group,'XTickLabel',[]);
newT = T;
for i=1:length(newList)
    newT(T==newList(i)) = i;
end
handles.T = newT;
choiceInput = sort(unique(T))';
set(handles.edit_cluster_list,'String',num2str(choiceInput));
linkaxes([handles.axes_posctrl  handles.axes_negctrl handles.axes_group handles.axes_selected],'x');
handles.binning = binning;
handles.newList = newList;
guidata(hObject, handles);  
handles = guidata(hObject); 

[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;
guidata(hObject, handles); 

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


% --- Executes on selection change in popupmenu_cluster.
function popupmenu_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_cluster contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_cluster
contents = cellstr(get(hObject,'String'));
switch get(hObject,'Value')
    case 1
        choiceInput = sort(unique(handles.T))';
        set(handles.edit_cluster_list,'String',num2str(choiceInput));
    otherwise
        cCluster = str2num(contents{get(hObject,'Value')});
        set(handles.edit_cluster_list,'String',num2str(cCluster));
end
[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;
guidata(hObject, handles); 


% --- Executes during object creation, after setting all properties.
function popupmenu_cluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_cluster (see GCBO)
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


% --- Executes on selection change in popupmenu_y_pca.
function popupmenu21_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_y_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_y_pca contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_y_pca


% --- Executes during object creation, after setting all properties.
function popupmenu21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_y_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plotparam.
function popupmenu_plotparam_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plotparam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plotparam contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plotparam


% --- Executes during object creation, after setting all properties.
function popupmenu_plotparam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plotparam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton_scatter3D.
function togglebutton_scatter3D_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_scatter3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit_commu,'String','3D Scatter plot');


% --- Executes on button press in togglebutton_contour.
function togglebutton_contour_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit_commu,'String','Contour plot');



% --- Executes on selection change in popupmenu23.
function popupmenu23_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu23 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu23


% --- Executes during object creation, after setting all properties.
function popupmenu23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_z_pca.
function popupmenu_z_pca_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_z_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_z_pca contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_z_pca


% --- Executes during object creation, after setting all properties.
function popupmenu_z_pca_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_z_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton_scatter2D.
function togglebutton_scatter2D_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_scatter2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit_commu,'String','2D Scatter plot');


function edit_selectedFields_Callback(hObject, eventdata, handles)
% hObject    handle to edit_selectedFields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_selectedFields as text
%        str2double(get(hObject,'String')) returns contents of edit_selectedFields as a double


% --- Executes during object creation, after setting all properties.
function edit_selectedFields_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_selectedFields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_selectedparams_Callback(hObject, eventdata, handles)
% hObject    handle to edit_selectedparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_selectedparams as text
%        str2double(get(hObject,'String')) returns contents of edit_selectedparams as a double
selectedparams = str2num(get(hObject,'String'));
table_data = get(handles.uitable_params,'Data');
param_names = [];

for i=1:size(table_data,1)
    param_names{i} = table_data{i,1};
end
Ind=1;
for i=selectedparams
    selectedparams_name{Ind} = param_names{i};
    Ind = Ind+1;
end
set(handles.popupmenu_selectedparams,'String',selectedparams_name);
set(handles.popupmenu_selectedparams,'Value',1);
guidata(hObject, handles);  

% --- Executes during object creation, after setting all properties.
function edit_selectedparams_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_selectedparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton_initialize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_initialize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in popupmenu_plotchoice.
function popupmenu_plotchoice_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plotchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plotchoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plotchoice
allInd = 1:size(handles.alldata,1);
switch get(hObject,'Value')
    case 1 % all
        handles.plotInd = allInd;
        handles.grayInd = [];
        set(handles.edit_phetype_list,'String',[]);
        guidata(hObject, handles);
        handles = guidata(hObject);
        [plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
        handles.plot_h = plot_h;
        handles.plotInd = plotInd;
        handles.grayInd = grayInd;
        
        
    case 2 % selected wells
        choiceInput = str2num(get(handles.edit_phetype_list,'String'));
        if isempty(choiceInput)
            choiceInput = sort(unique(handles.groupNo))';
            set(handles.edit_phetype_list,'String',num2str(choiceInput));
        end
        plotInd = [];
        for i=1:length(choiceInput)
            plotInd = [plotInd;find(handles.groupNo==choiceInput(i))];
        end
        handles.plotInd = sort(plotInd);
        handles.grayInd = sort(setdiff(allInd,plotInd));
    case 3 % selected phenotypes
        choiceInput = str2num(get(handles.edit_phetype_list,'String'));
        if isempty(choiceInput)
            choiceInput = sort(unique(handles.cellFate))';
            set(handles.edit_phetype_list,'String',num2str(choiceInput));
        end
        plotInd = [];
        for i=1:length(choiceInput)
            plotInd = [plotInd;find(handles.cellFate==choiceInput(i))];
        end
        handles.plotInd = sort(plotInd);
        handles.grayInd = sort(setdiff(allInd,plotInd));
    case 4 % selected cells
        if ~isempty(handles.selectedcellIndices)
            handles.plotInd = handles.selectedcellIndices;
            handles.grayInd = sort(setdiff(allInd,handles.selectedcellIndices));
        else
            set(handles.edit_commu,'String','No cell selected.');
        end
        
        set(handles.edit_phetype_list,'String',[]);
        guidata(hObject, handles);
        handles = guidata(hObject);
        [plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
        handles.plot_h = plot_h;
        handles.plotInd = plotInd;
        handles.grayInd = grayInd;
        
        
        
end
guidata(hObject, handles);
handles = guidata(hObject);
set(handles.edit_totalcellcount,'String',num2str(length(handles.plotInd)));



% --- Executes during object creation, after setting all properties.
function popupmenu_plotchoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plotchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_phetype_list_Callback(hObject, eventdata, handles)
% hObject    handle to edit_phetype_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_phetype_list as text
%        str2double(get(hObject,'String')) returns contents of edit_phetype_list as a double
[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;
guidata(hObject, handles); 


% --- Executes during object creation, after setting all properties.
function edit_phetype_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_phetype_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_totalcellcount_Callback(hObject, eventdata, handles)
% hObject    handle to edit_totalcellcount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_totalcellcount as text
%        str2double(get(hObject,'String')) returns contents of edit_totalcellcount as a double


% --- Executes during object creation, after setting all properties.
function edit_totalcellcount_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_totalcellcount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_plotIndividual.
function pushbutton_plotIndividual_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plotIndividual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

searchInd(1,:)= str2num(get(handles.edit_selectedRows,'String'));
searchInd(2,:)= str2num(get(handles.edit_selectedCols,'String'));
searchInd(3,:)= str2num(get(handles.edit_selectedFields,'String'));

outputsignalNo = str2num(get(handles.edit_outputsignalno,'String'));
dcm_obj = datacursormode(gcf);
infs = dcm_obj.getCursorInfo;

if ~isempty(infs)
    gind = getappdata(infs.Target,'gind');
    if ~isempty(gind)
        ind = infs.DataIndex;
        observationNo = gind(ind);
        
    elseif get(handles.popupmenu_c_pca,'Value') == 4
        c_X = get(handles.plot_h,'XData')';
        c_Y = get(handles.plot_h,'YData')';
        my_x = infs.Position(1);
        my_y = infs.Position(2);
        observationNo = find(c_X == my_x & c_Y == my_y);
        
        
    else
        set(handles.edit_commu,'String','Must select active (non-gray) point.');
        return;
    end
    
    if ~isempty(observationNo)
        currentData = handles.alldata(handles.plotInd,:);
        currentOriginalData = handles.originData(handles.plotInd,:);
        scell = currentData(observationNo,4);
        row = currentData(observationNo,1);
        col = currentData(observationNo,2);
        field = currentData(observationNo,3);
        H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
        peak_name  = ['/field' num2str(field)  '/peakmat' num2str(outputsignalNo)];
        fid = H5F.open(fullfile(handles.ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
        if H5L.exists(fid,peak_name,'H5P_DEFAULT')
            H5F.close(fid);
            p_peaks =  double(h5read(fullfile(handles.ndpathname,H5filename),peak_name,[double(scell) 1 1 1], [1 1 4 200]));
            t_peaks =  double(h5read(fullfile(handles.ndpathname,H5filename),peak_name,[double(scell) 2 1 1], [1 1 4 200]));
        end
        p_peaks = squeeze(p_peaks);
        p_truePeak = p_peaks(1,:);
        p_PeakHeight = p_peaks(2,:);
        p_PeakDuration = p_peaks(3,:);
        p_peakSelection = p_peaks(4,:);
        
        p_truePeak =p_truePeak(p_truePeak~=0);
        p_PeakHeight = p_PeakHeight(p_truePeak~=0);
        p_PeakDuration = p_PeakDuration(p_truePeak~=0);
        p_peakSelection = p_peakSelection(p_truePeak~=0);
        
        t_peaks = squeeze(t_peaks);
        t_truePeak = t_peaks(1,:);
        t_PeakHeight = t_peaks(2,:);
        t_PeakDuration = t_peaks(3,:);
        t_peakSelection = t_peaks(4,:);
        
        t_truePeak =t_truePeak(t_truePeak~=0);
        t_PeakHeight = t_PeakHeight(t_truePeak~=0);
        t_PeakDuration = t_PeakDuration(t_truePeak~=0);
        t_peakSelection = t_peakSelection(t_truePeak~=0);
        axes(handles.axes_individual);
        plot(handles.timestamp(currentOriginalData(observationNo,:)~=0),currentOriginalData(observationNo,currentOriginalData(observationNo,:)~=0),'k');
        YLim = get(handles.axes_individual,'YLim');
        for i=find(p_peakSelection==1)
            rectangle('Position',[p_truePeak(i),YLim(1),p_PeakDuration(i),YLim(2)-YLim(1)],...
                'FaceColor',[0.95 0.8 0.8],'EdgeColor','none','EraseMode','normal');hold on;
            
        end
        for i=find(t_peakSelection==1)
            rectangle('Position',[t_truePeak(i),YLim(1),t_PeakDuration(i),YLim(2)-YLim(1)],...
                'FaceColor',[ 0.8 0.8 0.95],'EdgeColor','none','EraseMode','normal');hold on;
        end
        plot(handles.timestamp(currentOriginalData(observationNo,:)~=0),currentOriginalData(observationNo,currentOriginalData(observationNo,:)~=0),'k');hold off;
        
        table_data = get(handles.uitable_params,'Data');
        for i=1:size(table_data,1)
            table_data{i,2} = nanmean(currentData(observationNo,i));
        end
        set(handles.uitable_params,'Data',table_data);
        set(handles.edit_commu,'String',['Showing cell# ' num2str(scell) ' from r' num2str(row) 'c' num2str(col) 'f' num2str(field)]);
    else
        set(handles.edit_commu,'String','Must select active (non-gray) point.');
    end
else
    set(handles.edit_commu,'String','Select a cell with cursor.');
end


% --------------------------------------------------------------------
function uitoggletool_showCurrentpoint_OffCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_zoomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uitoggletool_showCurrentpoint_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_zoomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit_commu,'String','Click data to show original signal trace. Only works when showing ''All'' data points.');


% --- Executes on button press in togglebutton_showselectedpolygon.
function togglebutton_showselectedpolygon_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_showselectedpolygon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_showselectedpolygon


if ~isempty(handles.originData) && ~isempty(handles.selectedcellIndices)
    old_xlim = get(handles.axes_pca,'XLim');
    old_ylim = get(handles.axes_pca,'YLim');
    [plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
    handles.plot_h = plot_h;
    handles.plotInd = plotInd;
    handles.grayInd = grayInd;
    set(handles.axes_pca,'XLim',old_xlim );
    set(handles.axes_pca,'YLim',old_ylim );
    if get(hObject,'Value') == 1
        searchInd(1,:)= str2num(get(handles.edit_selectedRows,'String'));
        searchInd(2,:)= str2num(get(handles.edit_selectedCols,'String'));
        searchInd(3,:)= str2num(get(handles.edit_selectedFields,'String'));
        mycolor = jet(size(searchInd,2));
        axes(handles.axes_individual);
        hold on;
        for i=1:length(handles.selectedcellIndices)
            plot(handles.timestamp(handles.originData(handles.selectedcellIndices(i),:)~=0),handles.originData(handles.selectedcellIndices(i),handles.originData(handles.selectedcellIndices(i),:)~=0),'Color',mycolor(handles.groupNo(handles.selectedcellIndices(i)),:));hold on;
        end
        hold off;
        drawnow;
    end
end




% --- Executes on button press in pushbutton_wellHistogram.
function pushbutton_wellHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_wellHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
searchInd(1,:)= str2num(get(handles.edit_selectedRows,'String'));
searchInd(2,:)= str2num(get(handles.edit_selectedCols,'String'));
searchInd(3,:)= str2num(get(handles.edit_selectedFields,'String'));
selectedInd = handles.selectedcellIndices;
wellPos = unique(handles.groupNo);
nelements = hist(handles.groupNo(selectedInd),wellPos);
bar(wellPos,nelements);
for i=1:length(wellPos)
    wellLabel{i} = ['r' num2str(searchInd(1,i)) 'c' num2str(searchInd(2,i))];
end
set(gca,'XTickLabel',wellLabel);



% --- Executes on button press in togglebutton_legendLogic.
function togglebutton_legendLogic_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_legendLogic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_legendLogic
switch get(hObject,'Value')
    case 1
        legend(handles.axes_pca,'show');
    case 0
        legend(handles.axes_pca,'hide');
end



function edit_nocontour_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nocontour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nocontour as text
%        str2double(get(hObject,'String')) returns contents of edit_nocontour as a double


% --- Executes during object creation, after setting all properties.
function edit_nocontour_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nocontour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_markersize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_markersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_markersize as text
%        str2double(get(hObject,'String')) returns contents of edit_markersize as a double


% --- Executes during object creation, after setting all properties.
function edit_markersize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_markersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_markertype_Callback(hObject, eventdata, handles)
% hObject    handle to edit_markertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_markertype as text
%        str2double(get(hObject,'String')) returns contents of edit_markertype as a double


% --- Executes during object creation, after setting all properties.
function edit_markertype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_markertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_contourcirclesize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_contourcirclesize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_contourcirclesize as text
%        str2double(get(hObject,'String')) returns contents of edit_contourcirclesize as a double


% --- Executes during object creation, after setting all properties.
function edit_contourcirclesize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_contourcirclesize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_meshsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_meshsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_meshsize as text
%        str2double(get(hObject,'String')) returns contents of edit_meshsize as a double


% --- Executes during object creation, after setting all properties.
function edit_meshsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_meshsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_outputsignalno_Callback(hObject, eventdata, handles)
% hObject    handle to edit_outputsignalno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_outputsignalno as text
%        str2double(get(hObject,'String')) returns contents of edit_outputsignalno as a double


% --- Executes during object creation, after setting all properties.
function edit_outputsignalno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_outputsignalno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton_orthorotate.
function togglebutton_orthorotate_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_orthorotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_orthorotate


% --- Executes on selection change in popupmenu_phenotype.
function popupmenu_phenotype_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_phenotype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_phenotype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_phenotype
contents = cellstr(get(hObject,'String'));
switch get(hObject,'Value')
    case 1
        choiceInput = sort(unique(handles.cellFate))';
        set(handles.edit_phetype_list,'String',num2str(choiceInput));
    otherwise
        regexp(contents{get(hObject,'Value')},'(?<=().expr(?=))');
        mycode = regexp(contents{get(hObject,'Value')},'.+\((?<code>.+)\)','tokens');
        set(handles.edit_phetype_list,'String',mycode{1}{1});
end
[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;
guidata(hObject, handles); 

% --- Executes during object creation, after setting all properties.
function popupmenu_phenotype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_phenotype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_wellposition.
function popupmenu_wellposition_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_wellposition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_wellposition contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_wellposition

searchInd(1,:)= str2num(get(handles.edit_selectedRows,'String'));
searchInd(2,:)= str2num(get(handles.edit_selectedCols,'String'));
searchInd(3,:)= str2num(get(handles.edit_selectedFields,'String'));
contents = cellstr(get(hObject,'String'));
switch get(hObject,'Value')
    case 1
        choiceInput = sort(unique(handles.groupNo))';
        set(handles.edit_position_list,'String',num2str(choiceInput));
        
    otherwise
        location = regexp(contents{get(hObject,'Value')},'r(?<row>\d+)c(?<col>\d+)','tokens');
        row = str2num(location{1}{1});
        col = str2num(location{1}{2});
        wells = sort(find(searchInd(1,:)==row & searchInd(2,:)==col));
        set(handles.edit_position_list,'String',num2str(wells));
end
[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;
guidata(hObject, handles); 

% --- Executes during object creation, after setting all properties.
function popupmenu_wellposition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_wellposition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_selection.
function popupmenu_selection_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_selection

[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;
guidata(hObject, handles); 

% --- Executes during object creation, after setting all properties.
function popupmenu_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton_reset.
function togglebutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_reset

set(handles.popupmenu_phenotype,'Value',1);
set(handles.popupmenu_wellposition,'Value',1);
set(handles.popupmenu_cluster,'Value',1);
set(handles.popupmenu_selection,'Value',1);
allInd = 1:size(handles.alldata,1);
handles.plotInd = allInd;
handles.grayInd = [];
guidata(hObject, handles);
handles = guidata(hObject);
[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;



function edit_position_list_Callback(hObject, eventdata, handles)
% hObject    handle to edit_position_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_position_list as text
%        str2double(get(hObject,'String')) returns contents of edit_position_list as a double
[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;
guidata(hObject, handles); 


% --- Executes during object creation, after setting all properties.
function edit_position_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_position_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cluster_list_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cluster_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cluster_list as text
%        str2double(get(hObject,'String')) returns contents of edit_cluster_list as a double
[plot_h,plotInd,grayInd] = plotPCA(handles.plottype,handles,handles.selectedcellIndices);
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;
guidata(hObject, handles); 


% --- Executes during object creation, after setting all properties.
function edit_cluster_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cluster_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in plotTypeChoice.
function plotTypeChoice_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in plotTypeChoice 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'togglebutton_scatter2D'
        [plot_h,plotInd,grayInd] = plotPCA(1,handles,handles.selectedcellIndices);
        handles.plottype = 1;
    case 'togglebutton_scatter3D'
        [plot_h,plotInd,grayInd] = plotPCA(2,handles,handles.selectedcellIndices);
        handles.plottype = 2;
    case 'togglebutton_contour'
        [plot_h,plotInd,grayInd] = plotPCA(3,handles,handles.selectedcellIndices);
        handles.plottype = 3;
end
handles.plot_h = plot_h;
handles.plotInd = plotInd;
handles.grayInd = grayInd;
guidata(hObject, handles);
