
function varargout = CellTracking(varargin)
%  varargin{1} - filetype 1 = PE export, 2 = tifstack, 3 = custom images
%  varargin{2} - trackinginfo  - [templateCH tp_1 tp_end]
%  varargin{3} - imagelocation - [row col field plane]
%  varargin{4} - channelnames eg. {'CFP';'FRET';'RFP'};
%  varargin{5} - fileformat or tiffstack file name
%  examples:
%  CellTracking(1,[templateCH tp_1 tp_end],[row col field plane])
%  CellTracking(2,[templateCH tp_1 tp_end],[],{'1';'2';'3'},<tiff stack file>)
%  CellTracking(3,[templateCH tp_1 tp_end],[],{'CFP';'YFP';'RFP'},'2012-11-16_%s_xy087_t%03g.tif')
%
% Last Modified by GUIDE v2.5 28-May-2013 23:07:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CellTracking_OpeningFcn, ...
    'gui_OutputFcn',  @CellTracking_OutputFcn, ...
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


function outputim = loadimage(filetype,fileformat,imlocation,tp,channelnames,SourceF)
row = imlocation(1);
col = imlocation(2);
field = imlocation(3);
plane = imlocation(4);
channel = imlocation(5);
totalCH = length(channelnames);
outputim = [];

switch filetype
    case 1
        
        filename = sprintf(fileformat,row,col,field,plane,channel,tp);
        if exist(fullfile(SourceF,filename),'file')
            outputim = imread(fullfile(SourceF,filename));
        end
    case 2
        if exist(fullfile(SourceF,fileformat),'file')
            outputim = imread(fullfile(SourceF,fileformat),'Index',totalCH*(tp-1)+channel);
        end
    case 3
        if ~isempty(channelnames)
            filename = sprintf(fileformat,channelnames{channel},tp);
        else
            filename = sprintf(fileformat,tp);
        end
        if exist(fullfile(SourceF,filename),'file');
            outputim = imread(fullfile(SourceF,filename));
        end
    case 4
        filename = sprintf(fileformat,row,col,field,tp,plane,channel);
        if exist(fullfile(SourceF,filename),'file')
            outputim = imread(fullfile(SourceF,filename));
        end  
end

% --- Executes just before CellTracking is made visible.
function CellTracking_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellTracking (see VARARGIN)

% Choose default command line output for CellTracking

handles.output = hObject;

filetype = [];
trackinginfo = [];
imagelocation = [];
channelnames = [];
fileformat = [];

switch length(varargin)
    case 1
        filetype = varargin{1};
    case 2
        filetype = varargin{1};
        if ischar(varargin{2}) && filetype == 3
            fileformat = varargin{2};
        else
            trackinginfo = varargin{2};  % [templateCH tp_1 tp_end]
        end
    case 3
        filetype = varargin{1};
        trackinginfo = varargin{2};
        imagelocation = varargin{3}; %[row col field plane]
    case 4
        filetype = varargin{1};
        trackinginfo = varargin{2};
        imagelocation = varargin{3};
        channelnames = varargin{4};  %{'CFP';'FRET';'RFP'};
    case 5
        filetype = varargin{1};
        trackinginfo = varargin{2};
        imagelocation = varargin{3};
        channelnames = varargin{4};
        fileformat = varargin{5};%
    case 0
        filetype = 0;
end
if isempty(filetype)
    display('Need to provide type (input#1) of input file(s)');
    return;
end



switch filetype
    case 1 % PE tiff input
        
        if isempty(imagelocation)
            display('Required image location (input#3) for PE-tiff format');
            return;
        end
        filetype = 1;
        set(handles.radiobutton_columbus,'Value',1);
        row = imagelocation(1);
        col = imagelocation(2);
        field = imagelocation(3);
        plane = imagelocation(4);
        
        if isempty(fileformat)
            fileformat = 'r%02.0fc%02.0ff%02.0fp%02.0frc%1.0f-ch1sk%ufk1fl1.tiff';
        end
        
        if ~isempty(channelnames)
            totalCHs = length(channelnames);
            set(handles.edit_totalCHs,'String',num2str(totalCHs));
        end
        
        channel = trackinginfo(1);
        set(handles.edit_CH,'String',num2str(channel));
        set(handles.edit_firstframe,'String',num2str(trackinginfo(2)));
        set(handles.edit_lastframe,'String',num2str(trackinginfo(3)));
        set(handles.edit_currentFrame,'String',num2str(trackinginfo(2)));
        
        
        handles.channelnames = channelnames;
        set(handles.edit_row,'String',num2str(row));
        set(handles.edit_col,'String',num2str(col));
        set(handles.edit_field,'String',num2str(field));
        set(handles.edit_plane,'String',num2str(plane));
        set(handles.edit_fileformat,'String',fileformat);
    case 2 % Tiff stack input
        
        if isempty(channelnames)
            display('Required channelnames (input#4) for tiffstack input');
            return;
        end
        
        if isempty(fileformat)
            display('Required filename (input#5) for tiffstack input');
            return;
        end
        filetype = 2;
        set(handles.radiobutton_tiffstack,'Value',1);
        
        totalCHs = length(channelnames);
        set(handles.edit_totalCHs,'String',num2str(totalCHs));
        [pathstr, name, ext] = fileparts(fileformat); %#ok<*NASGU,*ASGLU>
        if isempty(pathstr)
            fileformat = fullfile(pwd,fileformat);
        else
            fileformat = fullfile(pathstr,name);
        end
        
        if isempty(imagelocation)
            row = 1;
            col = 1;
            field = 1;
            plane = 1;
        else
            row = imagelocation(1);
            col = imagelocation(2);
            field = imagelocation(3);
            plane = imagelocation(4);
        end
        
        
        channel = trackinginfo(1);
        set(handles.edit_CH,'String',num2str(channel));
        set(handles.edit_firstframe,'String',num2str(trackinginfo(2)));
        set(handles.edit_lastframe,'String',num2str(trackinginfo(3)));
        set(handles.edit_currentFrame,'String',num2str(trackinginfo(2)));
        
        handles.channelnames = channelnames;
        set(handles.edit_row,'String',num2str(row));
        set(handles.edit_col,'String',num2str(col));
        set(handles.edit_field,'String',num2str(field));
        set(handles.edit_plane,'String',num2str(plane));
        set(handles.edit_fileformat,'String',fileformat);
    case 3 % tif file with custom naming input
        if isempty(fileformat)
            display('Required fileformat (input#5) for custom image input');
            return;
        end
        
        filetype = 3;
        set(handles.radiobutton_customtiff,'Value',1);
        
        if isempty(imagelocation)
            row = 1;
            col = 1;
            field = 1;
            plane = 1;
        else
            row = imagelocation(1);
            col = imagelocation(2);
            field = imagelocation(3);
            plane = imagelocation(4);
            
        end
        if strcmp(fileformat((end-2):end),'.nd')
            set(handles.edit_ndfilename,'String',fileformat);
            handles.ndfilename = fileformat;
            [notp stagePos stageName channelnames] = readndfile(fileformat);
            if notp==-1
                handles.initialframe = [];
                set(handles.edit_commu,'String',[fileformat ' does not exist. Please re-define input ND file.']);
                guidata(hObject, handles);
                return;
            end
            
            set(handles.edit_firstframe,'String',num2str(1));
            set(handles.edit_lastframe,'String',num2str(notp));
            set(handles.edit_currentFrame,'String',num2str(1));
            channel = 1;
            set(handles.popupmenu_stagePos,'String',stagePos);
            set(handles.popupmenu_stagePos,'Value',1);
            set(handles.edit_stageInfo,'String',stageName{1});
            handles.stageName = stageName;
            handles.channelnames = channelnames;
            
            prefix = fileformat(1:(end-3));
            handles.prefix = prefix;
            fileformat = [prefix '_%s_s1_t%g.tif'];
            set(handles.edit_fileformat,'String',fileformat);
            
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
        else
            set(handles.edit_fileformat,'String',fileformat);
            channel = trackinginfo(1);
            set(handles.edit_CH,'String',num2str(channel));
            set(handles.edit_firstframe,'String',num2str(trackinginfo(2)));
            set(handles.edit_lastframe,'String',num2str(trackinginfo(3)));
            set(handles.edit_currentFrame,'String',num2str(trackinginfo(2)));
        end
        
end

handles.filetype = filetype;
handles.cellpath = [];
handles.sisterList = [];
handles.res_cellpath = [];
handles.res_sisterList = [];
handles.SourceF = '/home/ss240/files/ImStor/sorger/data/NIC/Pat/';
guidata(hObject, handles);
handles = guidata(hObject);

cellsize = str2num(get(handles.edit_cellsize,'String')); %#ok<*ST2NM>
tp=str2num(get(handles.edit_currentFrame,'String'));
fftw('planner', 'hybrid');
if filetype~=0
    initialframe = loadimage(filetype,fileformat,[row col field plane channel],tp,channelnames,handles.SourceF);
    thres = stretchlim(initialframe);
    set(handles.edit_thresMin,'String',num2str(thres(1)));
    set(handles.edit_thresMax,'String',num2str(thres(2)));
    imshow(imadjust(initialframe,thres,[0 1]),'Parent',handles.axes1);
    maxF = str2num(get(handles.edit_lastframe,'String'));
    minF = str2num(get(handles.edit_firstframe,'String'));
    set(handles.slider_frames,'Max',maxF);
    set(handles.slider_frames,'Min',minF);
    set(handles.slider_frames,'Value',str2num(get(handles.edit_currentFrame,'String')));
    set(handles.slider_frames,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);
    set(handles.edit_commu,'String',{'Click AUTO for smart initialization. Click MANUAL for manual cell choosing (also possible after AUTO).'})
else
    initialframe = [];
    handles.filetype = 1;
    set(handles.edit_commu,'String',{'Please define input images to begin.'})
end

load fftexecutiontimes
handles.FFTiv=FFTiv;
handles.FFTrv=FFTrv;
handles.IFFTiv=IFFTiv;
handles.fcn1 = makeConstrainToRectFcn('impoint',get(handles.axes1,'XLim'),get(handles.axes1,'YLim'));

handles.cellpath = [];
handles.sisterList = [];
handles.initialframe = initialframe;
handles.greenflag=1;
handles.framestamp = 1;
handles.increment = 1;
handles.p = [];
handles.bg = [];
handles.bg_p = [];

guidata(hObject, handles);

% --- Executes on button press in pushbutton_auto.
function pushbutton_auto_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to pushbutton_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String','No input images.');
    return;
end

p = handles.p;
cellpath = handles.cellpath;
sisterList = handles.sisterList;
maxcellno = str2num(get(handles.edit_maxcells,'String'));

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
set(handles.togglebutton_editable,'Value',1);
tp = str2num(get(handles.edit_currentFrame,'String'));
currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
[coords_full] = autofindpoint(currentframe,handles.maxI,handles.invertedLog,str2num(get(handles.edit_cellareathres,'String')));
coords_max = sortrows([(1:size(coords_full,1))' (coords_full(:,1)-size(currentframe,2)/2).^2+(coords_full(:,2)-size(currentframe,1)/2).^2],2);
if maxcellno>size(coords_max,1)
    maxcellno = size(coords_max,1);
    set(handles.edit_maxcells,'String',num2str(maxcellno));
end

if get(handles.checkbox_randomized,'Value')
    RandInd = randperm(size(coords_max,1));
    coords  = coords_full(RandInd(1:maxcellno),:);
else
    coords  = coords_full(coords_max(1:maxcellno),:);
end
clear coords_full;
oldpsize = length(p);
for cell=1:size(coords,1)
    p{cell+oldpsize} = coordtopoint(coords(cell,:),cell,'r',handles.fcn1,handles.axes1,get(handles.checkbox_cellnostring,'Value'));
    currentPos = getPosition(p{cell});
    set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',coords(cell,1),coords(cell,2))));
    set(handles.edit_cellNo,'String',num2str(cell));
    addNewPositionCallback(p{cell+oldpsize},@(h) set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',h(1),h(2)))));
    addNewPositionCallback(p{cell+oldpsize},@(h) set(handles.edit_cellNo,'String',num2str(cell)));
    addNewPositionCallback(p{cell+oldpsize},@(h) set(handles.listbox_cells,'Value',cell));
    addNewPositionCallback(p{cell+oldpsize},@(h) set(handles.edit_sister,'String',num2str([-1 -1 -1])));
end

xy = zeros(length(p),2);
for i=1:length(p)
    xy(i,:) = getPosition(p{i});
end
cellpath{tp} = round([xy(:,1) xy(:,2)]);
sisterList{tp} = -1*ones(length(cellpath{tp}),3);

set(handles.listbox_cells,'String',num2str( [(1:size(cellpath{tp},1))' cellpath{tp} sisterList{tp}] ));
handles.cellpath = cellpath;
handles.sisterList = sisterList;
handles.p = p;
guidata(hObject, handles);
%set(handles.edit_commu,'String',['Recommended threshold = [' num2str(edgedThres(1)) ' ' num2str(edgedThres(2)) ']']);



% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
tp=str2num(get(handles.edit_currentFrame,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
maxF = str2num(get(handles.edit_lastframe,'String'));
minF = str2num(get(handles.edit_firstframe,'String'));
set(handles.slider_frames,'Max',maxF);
set(handles.slider_frames,'Min',minF);
set(handles.slider_frames,'Value',minF);
set(handles.edit_currentFrame,'String',num2str(minF));
set(handles.slider_frames,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);
if handles.filetype~=3
    handles.channelnames = [];
    guidata(hObject, handles);
    handles = guidata(hObject);
end
initialframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);

set(handles.edit_currentFrame,'String',num2str(tp));
set(handles.listbox_cells,'String',[]);
set(handles.listbox_bgs,'String',[]);
if ~isempty(initialframe)
    imshow(imadjust(initialframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
    set(handles.edit_commu,'String',{'Click AUTO for smart initialization. Click MANUAL for manual cell choosing (also possible after AUTO).'});
else
    set(handles.edit_commu,'String',{'Make sure the image path is correct.'});
end

handles.cellpath = [];
handles.bg = [];
handles.sisterList = [];
handles.res_cellpath = [];
handles.res_sisterList = [];
handles.initialframe = initialframe;
handles.greenflag=1;
handles.p = [];
handles.bg_p = [];
handles.framestamp = 1;
set(handles.radiobutton_frameno,'Value',1);
load fftexecutiontimes
handles.FFTiv=FFTiv;
handles.FFTrv=FFTrv;
handles.IFFTiv=IFFTiv;
handles.fcn1 = makeConstrainToRectFcn('impoint',get(handles.axes1,'XLim'),get(handles.axes1,'YLim'));
updateLists([],[],[],[],[],handles,tp);

guidata(hObject, handles);




% --- Executes on button press in pushbutton_deletecell.
function pushbutton_deletecell_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_deletecell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.initialframe)
    set(handles.edit_commu,'String','No input images.');
    return;
end

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));

cellpath = handles.cellpath;
bg = handles.bg;
sisterList = handles.sisterList;
cellpath = handles.cellpath;

selected_cell = str2num(get(handles.edit_cellNo,'String'));


c_tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

if ~isempty(cellpath)
    if length(cellpath)<last_tp
        last_tp = length(cellpath);
    end
end


for tp=first_tp:last_tp
    cellpath{tp}(selected_cell,:) = [-1 -1];
    sisterList{tp}(selected_cell,:) = [-1 -1 -1];
end
set(handles.edit_commu,'String',['Cell#' num2str(selected_cell) ' deleted']);

set(handles.togglebutton_editable,'Value',0);
currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],c_tp,handles.channelnames,handles.SourceF);
imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
[p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,c_tp,selected_cell);
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,c_tp);
handles.p = p;
handles.bg_p = bg_p;
handles.cellpath = cellpath;
handles.sisterList = sisterList;
guidata(hObject, handles);
drawnow;


% --- Executes on button press in pushbutton_manual.
function pushbutton_manual_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String','No input images.');
    return;
end

p = handles.p;
cellpath = handles.cellpath;
sisterList = handles.sisterList;
tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
set(handles.edit_commu,'String',{'Manually, Left-click to pick points and Right-click to choose the last point.'});
if isempty(cellpath)
    n=0;
else
    n=size(cellpath{tp},1);
end
but = 1;
while but == 1
    [xi,yi,but] = ginput(1);
    n = n+1;
    cellpath{tp}(n,:) = round([xi yi]);
    sisterList{tp}(n,:) = -1*ones(1,3);
    
    for t=first_tp:tp-1
        cellpath{t}(n,:) = [-1 -1];
        sisterList{t}(n,:) = -1*ones(1,3);
    end
    
    
    p{n} = coordtopoint([xi yi],n,'r',handles.fcn1,handles.axes1,get(handles.checkbox_cellnostring,'Value'));
    set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',xi,yi)));
    set(handles.edit_cellNo,'String',num2str(n));
    set(handles.listbox_cells,'String',num2str( [(1:size(cellpath{tp},1))' cellpath{tp} sisterList{tp}] ));
    set(handles.listbox_cells,'Value',n);
    addNewPositionCallback(p{n},@(h) set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',h(1),h(2)))));
    addNewPositionCallback(p{n},@(h) set(handles.edit_cellNo,'String',num2str(n)));
    addNewPositionCallback(p{n},@(h) set(handles.listbox_cells,'Value',n));
    addNewPositionCallback(p{n},@(h) set(handles.edit_sister,'String',num2str(sisterList{tp}(n,:))));
end

handles.cellpath = cellpath;
handles.sisterList = sisterList;
handles.p = p;
guidata(hObject, handles);

function p = coordtopoint(coord,id,color,fcn,axis,labellogic)
p = impoint(axis,coord);
setPositionConstraintFcn(p,fcn);
setColor(p,color);
if labellogic
    setString(p,num2str(id));
end


function [coords] = autofindpoint(M,maxI,invertLog,sizeThres)

%H = fspecial('unsharp');
%M = imfilter(M,H,'replicate');
if invertLog
    M = ((maxI)-1)-M;
end

minNucDiameter=sizeThres(1);
maxNucDiameter=sizeThres(2);
minFormfactor=0.5;
im=im2double(M);

maxNucArea=round(pi*maxNucDiameter^2/4);

combined_bw = zeros(size(M));
testSet = linspace(1.4,4.5,10);

for s = 1:length(testSet)
    
    [im2,thresh]  = edge(im,'canny',0,testSet(s));
    im3 = ~im2;
    bw=im3-bwareaopen(im3,maxNucArea,4);
    bw=imfill(bw,'holes');
    CC = bwconncomp(bw, 8);
    S = regionprops(CC,'EquivDiameter','Area','Perimeter');
    nucArea=cat(1,S.Area);
    nucPerim=cat(1,S.Perimeter);
    nucEquiDiameter=cat(1,S.EquivDiameter);
    nucFormfactor=4*pi*nucArea./(nucPerim.^2);
    L = labelmatrix(CC);
    BW2 = ismember(L, find(nucEquiDiameter >= minNucDiameter & ...
        nucEquiDiameter < maxNucDiameter & ...
        nucFormfactor>minFormfactor));
    
    %if sum(BW2(:)) > 0
    combined_bw = combined_bw | BW2;
    %end
end


%imshow(BW);
S = regionprops(combined_bw,'Centroid');
coords = zeros(length(S),2);
for i=1:length(S)
    coords(i,:) = S(i).Centroid;
end


function [x y BW] = templateToCentroid(M,xg,yg,maxI,invertLog)
if invertLog
    inverted = ((maxI)-1)-M;
end
BWc = zeros(size(M));
for i=1.2:0.6:3
    edgedIm = edge(M,'canny',0,i);
    BW = imfill(edgedIm,'holes');
    
    BW = bwmorph(BW,'open',1);
    BW = bwselect(BW,xg,yg);
    
    BWc = BWc | BW;
    
end

BW = BWc;
S  = regionprops(BW, 'centroid');

if isempty(find(BW==0, 1)) || isempty(find(BW==1, 1))
    x = xg;
    y = yg;
    
else
    x = round(S.Centroid(1));
    y = round(S.Centroid(2));
end

function [x y BW] = templateToCentroid2(M,xg,yg,maxI,invertLog)
if invertLog
    M = ((maxI)-1)-M;
end

BW = im2bw(imadjust(M), 1.3*graythresh(imadjust(M)));
BW = bwselect(BW,xg,yg);

S  = regionprops(BW, 'centroid');

if ~isfield(S,'Centroid') || length(S)>1
    x = xg;
    y = yg;
    
else
    x = round(S.Centroid(1));
    y = round(S.Centroid(2));
end


function [x y BW] = templateToCentroid3(M,xg,yg,maxI,invertLog,outerbox)
if invertLog
    M = ((maxI)-1)-M;
end
M = imadjust(M);
M(1,:) = 0;
M(end,:) = 0;
M(:,1) = 0;
M(:,end) = 0;
outBW = modchenvese(M,50,0.1,outerbox);
se = strel('disk',3);
outBW = imclose(outBW,se);
outBW(1,:) = 0;
outBW(end,:) = 0;
outBW(:,1) = 0;
outBW(:,end) = 0;

BW = bwselect(outBW,xg,yg);
S  = regionprops(BW, 'centroid');

if ~isfield(S,'Centroid') || length(S)>1
    x = xg;
    y = yg;
    
else
    x = round(S.Centroid(1));
    y = round(S.Centroid(2));
end

%   Adapted from code by Yue Wu (yue.wu@tufts.edu)
%   http://www.mathworks.com/matlabcentral/fileexchange/23445
function seg = modchenvese(I,num_iter,mu,outerbox)
I = imadjust(I);
s = outerbox./min(size(I,1),size(I,2)); % resize scale

n = zeros(size(I));
midY = round(size(n,1)/2);
midX = round(size(n,2)/2);
boxsize = round(size(I,1)*0.2);%round(outerbox/2*.9);
n(midY-boxsize:midY+boxsize,midX-boxsize:midX+boxsize) = 1;
I = imresize(I,s);
mask = imresize(n,s);

phi0 = bwdist(mask)-bwdist(1-mask)+im2double(mask)-.5;
force = eps;

%-- Display settings
%figure();
%subplot(2,2,1); imshow(I); title('Input Image');axis equal;
%subplot(2,2,2); contour(flipud(phi0), [0 0], 'r','LineWidth',1); title('initial contour');axis equal;
%subplot(2,2,3); title('Segmentation');axis equal;
%-- End Display original image and mask

%-- Main loop



for n=1:num_iter
    inidx = find(phi0>=0); % frontground index
    outidx = find(phi0<0); % background index
    force_image = 0; % initial image force for each layer
    L = im2double(I); % get one image component
    
    c1 = sum(sum(L.*Heaviside(phi0)))/(length(inidx)+eps); % average inside of Phi0
    c2 = sum(sum(L.*(1-Heaviside(phi0))))/(length(outidx)+eps); % verage outside of Phi0
    force_image=-(L-c1).^2+(L-c2).^2+force_image;
    % sum Image Force on all components (used for vector image)
    % if 'chan' is applied, this loop become one sigle code as a
    % result of layer = 1
    
    
    % calculate the external force of the image
    force = mu*kappa(phi0)./max(max(abs(kappa(phi0))))+force_image;
    
    % normalized the force
    force = force./(max(max(abs(force))));
    
    % get stepsize dt
    dt=0.5;
    
    % get parameters for checking whether to stop
    old = phi0;
    phi0 = phi0+dt.*force;
    new = phi0;
    indicator = checkstop(old,new,dt);
    
    
    
    % intermediate output
    %if(mod(n,20) == 0)
    %    showphi(I,phi0,n);
    %end;
    if indicator % decide to stop or continue
        %showphi(I,phi0,n);
        
        %make mask from SDF
        seg = phi0<=0; %-- Get mask from levelset
        
        midY = round(size(seg,1)/2);
        midX = round(size(seg,2)/2);
        if seg(midY,midX) > 0
            seg = bwselect(seg,midX,midY);
        else
            seg = ~seg;
            seg = bwselect(seg,midX,midY);
            %seg = bwselect(~seg,midX,midY);
        end
        
        
        %subplot(2,2,4); imshow(seg); title('Global Region-Based Segmentation');
        seg = imresize(seg,1/s);
        return;
    end
end;
%showphi(I,phi0,n);

%make mask from SDF
seg = phi0<=0; %-- Get mask from levelset

midY = round(size(seg,1)/2);
midX = round(size(seg,2)/2);

if seg(midY,midX) > 0
    seg = bwselect(seg,midX,midY);
else
    seg = ~seg;
    seg = bwselect(seg,midX,midY);
end


%subplot(2,2,4); imshow(seg); title('Global Region-Based Segmentation');
seg = imresize(seg,1/s);

function showphi(I, phi, i)
% show curve evolution of phi

% Copyright (c) 2009,
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved

for j = 1:size(phi,3)
    phi_{j} = phi(:,:,j); %#ok<*AGROW>
end
imshow(I,'initialmagnification','fit','displayrange',[0 255]);
hold on;

if size(phi,3) == 1
    contour(phi_{1}, [0 0], 'r','LineWidth',4);
    contour(phi_{1}, [0 0], 'g','LineWidth',1.3);
else
    contour(phi_{1}, [0 0], 'r','LineWidth',4);
    contour(phi_{1}, [0 0], 'x','LineWidth',1.3);
    contour(phi_{2}, [0 0], 'g','LineWidth',4);
    contour(phi_{2}, [0 0], 'x','LineWidth',1.3);
end
hold off;
title([num2str(i) ' Iterations']);
drawnow;

function KG = kappa(I)
% get curvature information of input image
% input: 2D image I
% output: curvature matrix KG

% Copyright (c) 2009,
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved

I = double(I);
[m,n] = size(I);
P = padarray(I,[1,1],1,'pre');
P = padarray(P,[1,1],1,'post');

% central difference
fy = P(3:end,2:n+1)-P(1:m,2:n+1);
fx = P(2:m+1,3:end)-P(2:m+1,1:n);
fyy = P(3:end,2:n+1)+P(1:m,2:n+1)-2*I;
fxx = P(2:m+1,3:end)+P(2:m+1,1:n)-2*I;
fxy = 0.25.*(P(3:end,3:end)-P(1:m,3:end)+P(3:end,1:n)-P(1:m,1:n));
G = (fx.^2+fy.^2).^(0.5);
K = (fxx.*fy.^2-2*fxy.*fx.*fy+fyy.*fx.^2)./((fx.^2+fy.^2+eps).^(1.5));
KG = K.*G;
KG(1,:) = eps;
KG(end,:) = eps;
KG(:,1) = eps;
KG(:,end) = eps;
KG = KG./max(max(abs(KG)));

function indicator = checkstop(old,new,dt)
% indicate whether we should performance further iteraions or stop

% Copyright (c) 2009,
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved

layer = size(new,3);

for i = 1:layer
    old_{i} = old(:,:,i);
    new_{i} = new(:,:,i);
end

if layer
    ind = find(abs(new)<=.5);
    M = length(ind);
    Q = sum(abs(new(ind)-old(ind)))./M;
    if Q<=dt*.18^2 %#ok<*BDSCI>
        indicator = 1;
    else
        indicator = 0;
    end
else
    ind1 = find(abs(old_{1})<1);
    ind2 = find(abs(old_{2})<1);
    M1 = length(ind1);
    M2 = length(ind2);
    Q1 = sum(abs(new_{1}(ind1)-old_{1}(ind1)))./M1;
    Q2 = sum(abs(new_{2}(ind2)-old_{2}(ind2)))./M2;
    if Q1<=dt*.18^2 && Q2<=dt*.18^2
        indicator = 1;
    else
        indicator = 0;
    end
end
return


function H=Heaviside(z)
% Heaviside step function (smoothed version)
% Copyright (c) 2009,
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved

Epsilon=10^(-5);
H=zeros(size(z,1),size(z,2));
idx1=find(z>Epsilon);
idx2=find(z<Epsilon & z>-Epsilon);
H(idx1)=1; %#ok<*FNDSB>
for i=1:length(idx2)
    H(idx2(i))=1/2*(1+z(idx2(i))/Epsilon+1/pi*sin(pi*z(idx2(i))/Epsilon));
end;



% --- Executes on button press in pushbutton_optimizeall.
function pushbutton_optimizeall_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_optimizeall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg=handles.bg;
cellsize = str2num(get(handles.edit_cellsize,'String'));
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));

fileformat = get(handles.edit_fileformat,'String');
tp = str2num(get(handles.edit_currentFrame,'String'));



currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
PosInd = find(cellpath{tp}(:,1)>0 & cellpath{tp}(:,2)>0)';
for cell=PosInd
    xL=max(cellpath{tp}(cell,1)-cellsize,1);
    xR=min(cellpath{tp}(cell,1)+cellsize,size(currentframe,2));
    yL=max(cellpath{tp}(cell,2)-cellsize,1);
    yR=min(cellpath{tp}(cell,2)+cellsize,size(currentframe,1));
    template = currentframe(yL:yR,xL:xR);
    [x1 y1 BW] = templateToCentroid(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog);
    cellpath{tp}(cell,:) = [xL+x1 yL+y1];
end
set(handles.edit_commu,'String',['Finished optimizing all points in frame:' num2str(tp)]);


%save points
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,tp);

handles.cellpath = cellpath;
guidata(hObject, handles);
handles = guidata(hObject);
restoreCellList(handles,tp);
imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;


% --- Executes on button press in pushbutton_optimize.
function pushbutton_optimize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_optimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.togglebutton_editable,'Value',1);
[p bg_p] = changeToEditable(handles);
handles.p = p;
handles.bg_p = bg_p;
guidata(hObject, handles);
handles = guidata(hObject);
cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg=handles.bg;
cellsize = str2num(get(handles.edit_cellsize,'String'));
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
tp = str2num(get(handles.edit_currentFrame,'String'));

currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);

cell = get(handles.listbox_cells,'Value');
xL=max(cellpath{tp}(cell,1)-cellsize,1);
xR=min(cellpath{tp}(cell,1)+cellsize,size(currentframe,2));
yL=max(cellpath{tp}(cell,2)-cellsize,1);
yR=min(cellpath{tp}(cell,2)+cellsize,size(currentframe,1));
template = currentframe(yL:yR,xL:xR);
[x1 y1 BW] = templateToCentroid(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog);
cellpath{tp}(cell,:) = [xL+x1 yL+y1];

%update the picture
clear template;
xL=max(cellpath{tp}(cell,1)-cellsize,1);
xR=min(cellpath{tp}(cell,1)+cellsize,size(currentframe,2));
yL=max(cellpath{tp}(cell,2)-cellsize,1);
yR=min(cellpath{tp}(cell,2)+cellsize,size(currentframe,1));
template = currentframe(yL:yR,xL:xR);
[x1 y1 BW] = templateToCentroid(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog);
edge = bwmorph(BW,'remove');
imshow(imadjust(template),'Parent',handles.axes2);
set(handles.axes2,'NextPlot','add');
[y x]=find(edge);
plot(handles.axes2,x,y,'.r','MarkerSize',1);
plot(handles.axes2,x1,y1,'xr');
set(handles.axes2,'NextPlot','replace');

%save new point
handles.cellpath = cellpath;
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,tp);
guidata(hObject, handles);
imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;

% --- Executes on button press in pushbutton_corrsetting.
function pushbutton_corrsetting_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_corrsetting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    cellpath = handles.cellpath;
    cellsize = str2num(get(handles.edit_cellsize,'String'));
    outersize = str2num(get(handles.edit_outersize,'String'));
    tp = str2num(get(handles.edit_currentFrame,'String'));
    cell = get(handles.listbox_cells,'Value');
    fcn = makeConstrainToRectFcn('imrect',get(handles.axes1,'XLim'),get(handles.axes1,'YLim'));
    h1 = imrect(handles.axes1, [cellpath{tp}(cell,1)-cellsize cellpath{tp}(cell,2)-cellsize 2*cellsize-1 2*cellsize-1]);
    
    setPositionConstraintFcn(h1,fcn);
    h2 = imrect(handles.axes1, [cellpath{tp}(cell,1)-outersize cellpath{tp}(cell,2)-outersize 2*outersize-1 2*outersize-1]);
    setResizable(h1,0);
    setResizable(h2,0);
    handles.rect1 = h1;
    handles.rect2 = h2;
elseif button_state == get(hObject,'Min')
    delete(handles.rect1);
    delete(handles.rect2);
end

guidata(hObject, handles);




% --- Executes on selection change in listbox_cells.
function listbox_cells_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_cells contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_cells
selected_cell = get(hObject,'Value');
handles.oldchosencell = selected_cell;
cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg=handles.bg;
c_tp = str2num(get(handles.edit_currentFrame,'String'));

first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));
if c_tp>last_tp
    tp=last_tp;
elseif c_tp<first_tp
    tp=first_tp;
else
    tp=c_tp;
end

if cellpath{tp}(selected_cell,1) ~= -1 || cellpath{tp}(selected_cell,2) ~= -1
    
    
    set(handles.edit_cellNo,'String',num2str(selected_cell));
    set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',cellpath{tp}(selected_cell,1),cellpath{tp}(selected_cell,2))));
    set(handles.edit_sister,'String',num2str(sisterList{tp}(selected_cell,:)));
    
    cellsize = str2num(get(handles.edit_cellsize,'String'));
    row = str2num(get(handles.edit_row,'String'));
    col = str2num(get(handles.edit_col,'String'));
    field = str2num(get(handles.edit_field,'String'));
    plane = str2num(get(handles.edit_plane,'String'));
    channel= str2num(get(handles.edit_CH,'String'));
    
    
    currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
    set(handles.edit_currentFrame,'String',num2str(tp));
    imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
    if get(handles.checkbox_cellmarking,'Value')
        [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
        handles.p = p;
        handles.bg_p = bg_p;
    end
    
    xL=max(cellpath{tp}(selected_cell,1)-cellsize,1);
    xR=min(cellpath{tp}(selected_cell,1)+cellsize,size(currentframe,2));
    yL=max(cellpath{tp}(selected_cell,2)-cellsize,1);
    yR=min(cellpath{tp}(selected_cell,2)+cellsize,size(currentframe,1));
    template = currentframe(yL:yR,xL:xR);
    
    switch handles.cellrecogmethod
        case 1
            [x1 y1 BW] = templateToCentroid(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog);
        case 2
            [x1 y1 BW] = templateToCentroid2(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog);
        case 3
            [x1 y1 BW] = templateToCentroid3(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog,str2num(get(handles.edit_outersize,'String')));

    end
    
    edge = bwmorph(BW,'remove');
    imshow(imadjust(template),'Parent',handles.axes2);
    set(handles.axes2,'NextPlot','add');
    [y x]=find(edge);
    plot(handles.axes2,x,y,'.r','MarkerSize',1);
    plot(handles.axes2,cellpath{tp}(selected_cell,1)-xL,cellpath{tp}(selected_cell,2)-yL,'xr');
    set(handles.axes2,'NextPlot','replace');
    guidata(hObject, handles);
    drawnow;
    set(handles.edit_commu,'String',['Showing Cell#' num2str(selected_cell)]);
else
    set(handles.edit_cellNo,'String',num2str(selected_cell));
    set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',cellpath{tp}(selected_cell,1),cellpath{tp}(selected_cell,2))));
    set(handles.edit_sister,'String',num2str(sisterList{tp}(selected_cell,:)));
    
    set(handles.edit_commu,'String',['Cell#' num2str(selected_cell) ' does not exist at current time point.']);
end
% --- Outputs from this function are returned to the command line.
function varargout = CellTracking_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_recordpositions.
function pushbutton_recordpositions_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_recordpositions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.togglebutton_editable,'Value') == get(handles.togglebutton_editable,'Max')
    p=handles.p;
    cellpath = handles.cellpath;
    sisterList = handles.sisterList;
    bg = handles.bg;
    bg_p = handles.bg_p;
    selected_cell = str2num(get(handles.edit_cellNo,'String'));
    c_tp = str2num(get(handles.edit_currentFrame,'String'));
    
    for i=1:length(p)
        if ~isempty(p{i})
            new_xy = getPosition(p{i});
            if ~isempty(cellpath) && ~isempty(cellpath{c_tp})
                cellpath{c_tp}(i,:) = round(new_xy);
            end
        end
    end
    
    for i=1:length(bg_p)
        xy_bg(i,:) = getPosition(bg_p{i});
    end
    
    if ~isempty(bg) && ~isempty(bg{c_tp})
        bg{c_tp} = round([xy_bg(:,1) xy_bg(:,2)]);
        set(handles.listbox_bgs,'String',num2str( [(1:size(bg{c_tp},1))' bg{c_tp}] ));
    end
    
    set(handles.edit_commu,'String',['Recorded data to frame ' num2str(c_tp) '.']);
    handles.cellpath = cellpath;
    handles.bg = bg;
    updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,c_tp);
    guidata(hObject, handles);
end
% --- Executes on button press in pushbutton_tofirstframe.
function pushbutton_tofirstframe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tofirstframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String','No input images.');
    return;
end

set(handles.edit_commu,'String',{'Initial frame'});
cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg = handles.bg;

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
fileformat = get(handles.edit_fileformat,'String');
tp = str2num(get(handles.edit_firstframe,'String'));
set(handles.slider_frames,'Value',tp);
firstframe = loadimage(handles.filetype,fileformat,[row col field plane channel],tp,handles.channelnames,handles.SourceF);
set(handles.edit_currentFrame,'String',num2str(tp));

set(handles.slider_frames,'Value',tp);
if handles.filetype == 3
    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},first_tp);
    else
        filename = sprintf(fileformat,first_tp);
    end
    first_info = imfinfo(fullfile(handles.SourceF,filename));
    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},tp);
    else
        filename = sprintf(fileformat,tp);
    end
    current_info = imfinfo(fullfile(handles.SourceF,filename));
    [Y, M, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
    hour = 24*D+round(H);
    minute = round(MN);
    second = round(S);
    set(handles.edit_commu,'String',[num2str(hour,'%02.0f') ':' num2str(minute,'%02.0f') ':' num2str(second,'%02.0f')]);
else
    set(handles.edit_commu,'String',{'First frame'});
end

imshow(imadjust((firstframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,tp);

function [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,res_cellpath,res_sisterList,bg,tp,selected_cell)
p = [];
bg_p = [];

if get(handles.togglebutton_editable,'Value')
    
    if ~isempty(cellpath) && ~isempty(cellpath{tp})
        PosInd = find(cellpath{tp}(:,1)>0)';
        for cell=PosInd
            if cell == selected_cell
                p{cell} = coordtopoint(cellpath{tp}(cell,:),cell,'b',handles.fcn1,handles.axes1,get(handles.checkbox_cellnostring,'Value'));
                addNewPositionCallback(p{cell},@(h) set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',h(1),h(2)))));
                addNewPositionCallback(p{cell},@(h) set(handles.edit_cellNo,'String',num2str(cell)));
                addNewPositionCallback(p{cell},@(h) set(handles.listbox_cells,'Value',cell));
                addNewPositionCallback(p{cell},@(h) set(handles.edit_sister,'String',num2str(sisterList{tp}(cell,:))));
            end
        end
        
    end
    
    if ~isempty(bg)
        for b=1:size(bg{tp},1)
            bg_p{b} = coordtopoint(bg{tp}(b,:),b,'g',handles.fcn1,handles.axes1,get(handles.checkbox_cellnostring,'Value'));
            addNewPositionCallback(bg_p{b},@(h) set(handles.listbox_bgs,'Value',b));
        end
        
    end
    
else
    hold on;
    axes(handles.axes1); hold on;
    set(handles.edit_currentFrame,'String',num2str(tp));
    if ~isempty(cellpath) && ~isempty(cellpath{tp})
        PosInd = find(cellpath{tp}(:,1)>0);
        nosisterInd = find(sisterList{tp}(PosInd,1)==-1);
        sister1Ind = find(sisterList{tp}(PosInd,1)~=-1 & sisterList{tp}(PosInd,2)==-1 & sisterList{tp}(PosInd,3)==-1 );
        sister2Ind = find(sisterList{tp}(PosInd,1)~=-1 & sisterList{tp}(PosInd,2)~=-1 & sisterList{tp}(PosInd,3)==-1 );
        sister3Ind = find(sisterList{tp}(PosInd,1)~=-1 & sisterList{tp}(PosInd,2)~=-1 & sisterList{tp}(PosInd,3)~=-1 );
        
        plot(cellpath{tp}(PosInd(nosisterInd),1),cellpath{tp}(PosInd(nosisterInd),2),'x','MarkerFaceColor','c','MarkerEdgeColor','c'); hold on;
        plot(cellpath{tp}(PosInd(sister1Ind),1),cellpath{tp}(PosInd(sister1Ind),2),'o','MarkerFaceColor','r','MarkerEdgeColor','r'); hold on;
        plot(cellpath{tp}(PosInd(sister2Ind),1),cellpath{tp}(PosInd(sister2Ind),2),'o','MarkerFaceColor','g','MarkerEdgeColor','g'); hold on;
        plot(cellpath{tp}(PosInd(sister3Ind),1),cellpath{tp}(PosInd(sister3Ind),2),'o','MarkerFaceColor','b','MarkerEdgeColor','b'); hold on;
        
        if get(handles.checkbox_cellnostring,'Value')
            text(cellpath{tp}(PosInd,1)+5,cellpath{tp}(PosInd,2)+5,num2str(PosInd),'HorizontalAlignment','left',...
                'VerticalAlignment','middle','color',[0 .9 .5]);
        end
        if selected_cell>size(cellpath{tp},1)
            selected_cell = 1;
        end
        plot(cellpath{tp}(selected_cell,1),cellpath{tp}(selected_cell,2),'o','MarkerFaceColor','none','MarkerEdgeColor','y'); hold on;
        if get(handles.checkbox_cellnostring,'Value') ~=1
            text(cellpath{tp}(selected_cell,1)+10,cellpath{tp}(selected_cell,2)+5,num2str(selected_cell),'HorizontalAlignment','left',...
                'VerticalAlignment','middle','color',[1 1 0]);
        end
        
    end
    if ~isempty(bg)
        plot(bg{tp}(:,1),bg{tp}(:,2),'.','MarkerFaceColor','c','MarkerEdgeColor','c');
    end
    hold off;
end



if ~isempty(res_cellpath) && ~isempty(res_cellpath{tp})
    hold on;
    plot(res_cellpath{tp}(:,1),res_cellpath{tp}(:,2),'x','MarkerFaceColor','w','MarkerEdgeColor','w'); hold on;
    if get(handles.checkbox_cellnostring,'Value')
        text(res_cellpath{tp}(:,1)+5,res_cellpath{tp}(:,2)+5,num2str((1:length(res_cellpath{tp}))'),'HorizontalAlignment','left',...
            'VerticalAlignment','middle','color','w');
    end
    if isempty(cellpath) || isempty(cellpath{tp})
        if selected_cell>size(res_cellpath{tp},1)
            selected_cell = 1;
        end
        plot(res_cellpath{tp}(selected_cell,1),res_cellpath{tp}(selected_cell,2),'o','MarkerFaceColor','none','MarkerEdgeColor','y'); hold on;
        if get(handles.checkbox_cellnostring,'Value') ~=1
            text(res_cellpath{tp}(selected_cell,1)+10,res_cellpath{tp}(selected_cell,2)+5,num2str(selected_cell),'HorizontalAlignment','left',...
                'VerticalAlignment','middle','color',[1 1 0]);
        end
    end
    hold off;
end


function updateLists(cellpath,sisterList,res_cellpath,res_sisterList,bg,handles,tp)

if ~isempty(cellpath) && ~isempty(cellpath{tp})
    set(handles.listbox_cells,'String',num2str( [(1:size(cellpath{tp},1))' cellpath{tp} sisterList{tp}] ));
else
    set(handles.listbox_cells,'String',[]);
end

if ~isempty(res_cellpath) && ~isempty(res_cellpath{tp})
    set(handles.listbox_Restcells,'String',num2str( [(1:size(res_cellpath{tp},1))' res_cellpath{tp} res_sisterList{tp}] ));
else
    set(handles.listbox_Restcells,'String',[]);
end
if ~isempty(bg) && ~isempty(bg{tp})
    set(handles.listbox_bgs,'String',num2str( [(1:size(bg{tp},1))' bg{tp}] ));
else
    set(handles.listbox_bgs,'String',[]);
end

% --- Executes on button press in pushbutton_topreviousframe.
function pushbutton_topreviousframe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_topreviousframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String','No input images.');
    return;
end
cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg=handles.bg;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
c_tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
fileformat = get(handles.edit_fileformat,'String');
if c_tp-1<first_tp
    pre_tp=first_tp;
else
    pre_tp=c_tp-1;
end

set(handles.slider_frames,'Value',pre_tp);
previousframe = loadimage(handles.filetype,fileformat,[row col field plane channel],pre_tp,handles.channelnames,handles.SourceF);
if handles.filetype == 3
    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},first_tp);
    else
        filename = sprintf(fileformat,first_tp);
    end
    first_info = imfinfo(fullfile(handles.SourceF,filename));
    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},pre_tp);
    else
        filename = sprintf(fileformat,pre_tp);
    end
    current_info = imfinfo(fullfile(handles.SourceF,filename));
    [Y, M, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
    hour = 24*D+round(H);
    minute = round(MN);
    second = round(S);
    set(handles.edit_commu,'String',[num2str(hour,'%02.0f') ':' num2str(minute,'%02.0f') ':' num2str(second,'%02.0f')]);
end

set(handles.edit_currentFrame,'String',num2str(pre_tp));
imshow(imadjust((previousframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,pre_tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,pre_tp);
% --- Executes on button press in pushbutton_tonextframe.
function pushbutton_tonextframe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tonextframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String','No input images.');
    return;
end

cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg=handles.bg;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
c_tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));
fileformat = get(handles.edit_fileformat,'String');
if c_tp+1>last_tp
    post_tp=last_tp;
else
    post_tp=c_tp+1;
end
if handles.filetype == 3
    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},first_tp);
    else
        filename = sprintf(fileformat,first_tp);
    end
    first_info = imfinfo(fullfile(handles.SourceF,filename));
end
set(handles.slider_frames,'Value',post_tp);
nextframe = loadimage(handles.filetype,fileformat,[row col field plane channel],post_tp,handles.channelnames,handles.SourceF);
if handles.filetype == 3
    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},post_tp);
    else
        filename = sprintf(fileformat,post_tp);
    end
    current_info = imfinfo(fullfile(handles.SourceF,filename));
    [Y, M, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
    hour = 24*D+round(H);
    minute = round(MN);
    second = round(S);
    set(handles.edit_commu,'String',[num2str(hour,'%02.0f') ':' num2str(minute,'%02.0f') ':' num2str(second,'%02.0f')]);
end
imshow(imadjust((nextframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
set(handles.edit_currentFrame,'String',num2str(post_tp));

if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,post_tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,post_tp);
% --- Executes on button press in pushbutton_tolastframe.
function pushbutton_tolastframe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tolastframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String','No input images.');
    return;
end


cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg=handles.bg;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
tp = str2num(get(handles.edit_lastframe,'String'));
set(handles.slider_frames,'Value',tp);
first_tp = str2num(get(handles.edit_firstframe,'String'));
fileformat = get(handles.edit_fileformat,'String');
lastframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
imshow(imadjust((lastframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
set(handles.edit_currentFrame,'String',num2str(tp));


if handles.filetype == 3

    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},first_tp);
    else
        filename = sprintf(fileformat,first_tp);
    end
    first_info = imfinfo(fullfile(handles.SourceF,filename));
    
    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},tp);
    else
        filename = sprintf(fileformat,tp);
    end
    current_info = imfinfo(fullfile(handles.SourceF,filename));
    [Y, M, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
    hour = 24*D+round(H);
    minute = round(MN);
    second = round(S);
    set(handles.edit_commu,'String',[num2str(hour,'%02.0f') ':' num2str(minute,'%02.0f') ':' num2str(second,'%02.0f')]);
else
    set(handles.edit_commu,'String',{'Final frame'});
end


if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,tp);

function edit_currentFrame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_currentFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_currentFrame as text
%        str2double(get(hObject,'String')) returns contents of edit_currentFrame as a double
if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']); %#ok<*NBRAK>
    return;
end

cellpath = handles.cellpath;
bg = handles.bg;
sisterList = handles.sisterList;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
c_tp = str2num(get(hObject,'String'));
set(handles.slider_frames,'Value',str2num(get(hObject,'String')));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));
if c_tp>last_tp
    tp=last_tp;
elseif c_tp<first_tp
    tp=first_tp;
else
    tp=c_tp;
end
set(handles.edit_currentFrame,'String',num2str(tp));

currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,tp);

function restoreCellList(handles,tp)
cellpath = handles.cellpath;
sisterList = handles.sisterList;

selected_cell = get(handles.listbox_cells,'Value');
set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',cellpath{tp}(selected_cell,1),cellpath{tp}(selected_cell,2))));
set(handles.edit_cellNo,'String',num2str(selected_cell));
set(handles.edit_sister,'String',num2str(sisterList{tp}(selected_cell,:)));


% --- Executes on button press in pushbutton_play.
function pushbutton_play_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg=handles.bg;

fileformat = get(handles.edit_fileformat,'String');
c_tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

switch handles.increment
    case 1
        endFrame = last_tp;
        for tp=first_tp:1:last_tp
            if ~isempty(cellpath)
                if length(cellpath)>=tp
                    endFrame = tp;
                else
                    break;
                end
                
            end
        end
    case -1
        endFrame   = first_tp;
        for tp=last_tp:-1:first_tp
            if ~isempty(cellpath) && ~isempty(cellpath{tp})
                endFrame = tp;
            else
                break;
            end
        end
end



if c_tp~=endFrame+handles.increment
    set(handles.edit_commu,'String',{'Playing...'});
    handles.greenflag = 1;
    guidata(hObject, handles);
    handles = guidata(hObject);
    
    if handles.filetype == 3
        if ~isempty(handles.channelnames)
            filename = sprintf(fileformat,handles.channelnames{channel},first_tp);
        else
            filename = sprintf(fileformat,first_tp);
        end

        first_info = imfinfo(fullfile(handles.SourceF,filename));
    end
    for tp=c_tp:handles.increment:endFrame
        set(handles.slider_frames,'Value',tp);
        currentframe = loadimage(handles.filetype,fileformat,[row col field plane channel],tp,handles.channelnames,handles.SourceF);
        if handles.filetype == 3
            if ~isempty(handles.channelnames)
                filename = sprintf(fileformat,handles.channelnames{channel},tp);
            else
                filename = sprintf(fileformat,tp);
            end
            current_info = imfinfo(fullfile(handles.SourceF,filename));
            [Y, M, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
            hour = 24*D+round(H);
            minute = round(MN);
            second = round(S);
            set(handles.edit_commu,'String',[num2str(hour,'%02.0f') ':' num2str(minute,'%02.0f') ':' num2str(second,'%02.0f')]);
        end
        
        imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
        set(handles.edit_currentFrame,'String',num2str(tp));
        
        if get(handles.checkbox_cellmarking,'Value')
            [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
            handles.p = p;
            handles.bg_p = bg_p;
            guidata(hObject, handles);
        end
        drawnow;
        updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,tp);
        handles = guidata(hObject);
        if handles.greenflag==0
            break;
        end
        
    end
    if ~isempty(cellpath) && ~isempty(cellpath{tp})
        set(handles.listbox_cells,'String',num2str( [(1:size(cellpath{tp},1))' cellpath{tp} sisterList{tp}] ));
    end
    handles.greenflag = 1;
    guidata(hObject, handles);
end


% --- Executes on button press in pushbutton_pause.
function pushbutton_pause_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end

handles.greenflag = 0;
guidata(hObject, handles);
set(handles.edit_commu,'String',{'Paused'});


% --- Executes on button press in togglebutton_forward.
function togglebutton_forward_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg = handles.bg;

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));

currentFrame = str2num(get(handles.edit_currentFrame,'String'));
lastframe   = str2num(get(handles.edit_lastframe,'String'));

nextframe = currentFrame+1;
if nextframe>lastframe
    nextframe=lastframe;
end
if nextframe == currentFrame
    return;
end


FFTiv=handles.FFTiv;
FFTrv=handles.FFTrv;
IFFTiv=handles.IFFTiv;


previousframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],currentFrame,handles.channelnames,handles.SourceF);
currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],nextframe,handles.channelnames,handles.SourceF);

opt = detbestlength2(FFTrv,FFTiv,IFFTiv,size(imresize(currentframe,0.1)),size(imresize(previousframe,0.1)),isreal(imresize(currentframe,0.1)),isreal(imresize(previousframe,0.1)));
[x y maxVal] = corrMatching2(imresize(currentframe,0.1),imresize(previousframe,0.1),opt);

xshift = 10*x-round(size(currentframe,2)/2);
yshift = 10*y-round(size(currentframe,1)/2);

cellpath{nextframe} = [cellpath{currentFrame}(:,1)+xshift cellpath{currentFrame}(:,2)+yshift];
sisterList{nextframe} = sisterList{currentFrame};

if ~isempty(bg)
    bg{nextframe} = [bg{currentFrame}(:,1)+xshift bg{currentFrame}(:,2)+yshift];
    set(handles.listbox_bgs,'String',num2str( [(1:size(bg{nextframe},1))' bg{nextframe}] ));
end
set(handles.listbox_cells,'String',num2str( [(1:size(cellpath{nextframe},1))' cellpath{nextframe} sisterList{nextframe}] ));
handles.cellpath = cellpath;
handles.sisterList = sisterList;
handles.bg = bg;
guidata(hObject, handles);
handles = guidata(hObject);
restoreCellList(handles,nextframe);
imshow(imadjust((currentframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
set(handles.edit_currentFrame,'String',num2str(nextframe));
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,nextframe,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;

% --- Executes when selected object is changed in uibuttongroup_cellrecogmethod.
function uibuttongroup_cellrecogmethod_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup_cellrecogmethod
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton_method1'
        set(handles.edit_commu,'String','Canny Edging chosen');
        handles.cellrecogmethod = 1;
    case 'radiobutton_method2'
        set(handles.edit_commu,'String','Otsu Thresholding');
        handles.cellrecogmethod = 2;
    case 'radiobutton_method3'
        set(handles.edit_commu,'String','Contour tracking');
        handles.cellrecogmethod = 3;
    otherwise
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton_correlationthru.
function pushbutton_correlationthru_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_correlationthru (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end
set(handles.togglebutton_editable,'Value',0);
cellpath = handles.cellpath;
sisterList = handles.sisterList;
res_cellpath = handles.res_cellpath;
res_sisterList = handles.res_sisterList;
bg=handles.bg;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
outersize = str2num(get(handles.edit_outersize,'String'));

FFTiv=handles.FFTiv;
FFTrv=handles.FFTrv;
IFFTiv=handles.IFFTiv;
cFrame = str2num(get(handles.edit_currentFrame,'String'));
switch handles.increment
    case 1
        endFrame   = str2num(get(handles.edit_lastframe,'String'))-1;
    case -1
        endFrame   = str2num(get(handles.edit_firstframe,'String'))+1;
end

if cFrame~=endFrame+handles.increment
    
    opt = detbestlength2(FFTrv,FFTiv,IFFTiv,2*[outersize outersize],2*[cellsize cellsize],1,1);
    
    handles.greenflag = 1;
    guidata(hObject, handles);
    handles = guidata(hObject);
    
    set(handles.edit_commu,'String',{'Processing...'});
    for t=cFrame:handles.increment:endFrame
        
        tp=t;
        previousframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
        tp=t+handles.increment;
        currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
        imshow(imadjust((currentframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
        
        
        set(handles.edit_currentFrame,'String',num2str(tp));
        set(handles.slider_frames,'Value',tp);
        
        if length(sisterList) >= tp && ~isempty(sisterList{tp}) && ~isempty(sisterList{1})
            sisExistInd = find(sisterList{tp}(:,1) ~= -1 & sisterList{tp}(:,1) ~= 0);
        else
            sisExistInd = [];
        end
        
        PosInd = find(cellpath{tp-handles.increment}(:,1)>0 & cellpath{tp-handles.increment}(:,2)>0)';
        
        for cell=PosInd
            if isempty(find(cell==sisExistInd,1))
                sisterList{tp}(cell,:) = sisterList{tp-handles.increment}(cell,:);
            end
            
            xL=max(cellpath{tp-handles.increment}(cell,1)-cellsize,1);
            xR=min(cellpath{tp-handles.increment}(cell,1)+cellsize-1,size(previousframe,2));
            yL=max(cellpath{tp-handles.increment}(cell,2)-cellsize,1);
            yR=min(cellpath{tp-handles.increment}(cell,2)+cellsize-1,size(previousframe,1));
            
            template = previousframe(yL:yR,xL:xR);
            
            xL=max(cellpath{tp-handles.increment}(cell,1)-outersize,1);
            xR=min(cellpath{tp-handles.increment}(cell,1)+outersize-1,size(currentframe,2));
            yL=max(cellpath{tp-handles.increment}(cell,2)-outersize,1);
            yR=min(cellpath{tp-handles.increment}(cell,2)+outersize-1,size(currentframe,1));
            
            testframe = currentframe(yL:yR,xL:xR);
            
            [x1 y1 maxVal] = corrMatching2(testframe, template,opt);
            if isempty(x1) || maxVal<str2num(get(handles.edit_similarityThres,'String'))
                x1 = cellpath{tp-handles.increment}(cell,1)-xL;
                y1 = cellpath{tp-handles.increment}(cell,2)-yL;
            end
            
            if get(handles.checkbox_autooptimize,'Value')==1
                [x1 y1 BW] = templateToCentroid(testframe,x1,y1,handles.maxI,handles.invertedLog);
            end
            
            cellpath{tp}(cell,:) = [xL+x1 yL+y1];
        end
        
        DeathInd = find(cellpath{tp-handles.increment}(:,1)==-2)';
        
        for cell=DeathInd
            sisterList{tp}(cell,:) = sisterList{tp-handles.increment}(cell,:);
            cellpath{tp}(cell,:) = cellpath{tp-handles.increment}(cell,:);
        end
        
        if ~isempty(bg)
            %updating bg points by
            xshift = mean(cellpath{tp}(:,1))-mean(cellpath{tp-handles.increment}(:,1));
            yshift = mean(cellpath{tp}(:,2))-mean(cellpath{tp-handles.increment}(:,2));
            bg{tp}(:,1) = bg{tp-handles.increment}(:,1)+round(xshift);
            bg{tp}(:,2) = bg{tp-handles.increment}(:,2)+round(yshift);
            handles.bg = bg;
        end
        updateLists(cellpath,sisterList,res_cellpath,res_sisterList,bg,handles,tp);
        guidata(hObject, handles);
        handles = guidata(hObject);
        if get(handles.checkbox_cellmarking,'Value')
            [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
            handles.p = p;
            handles.bg_p = bg_p;
            guidata(hObject, handles);
        end
        drawnow;
        handles = guidata(hObject);
        if handles.greenflag==0
            break;
        end
    end
    handles.greenflag = 1;
    handles.cellpath = cellpath;
    handles.sisterList = sisterList;
    if tp==endFrame
        set(handles.edit_commu,'String',{'Optimize points and press Record to save positions for each frame.'});
    end
    guidata(hObject, handles);
    handles = guidata(hObject);
    restoreCellList(handles,tp);
end
% --- Executes on button press in pushbutton_correlation.
function pushbutton_correlation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_correlation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end
set(handles.togglebutton_editable,'Value',0);
cellpath = handles.cellpath;
sisterList = handles.sisterList;
res_cellpath = handles.res_cellpath;
res_sisterList = handles.res_sisterList;
bg = handles.bg;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
outersize = str2num(get(handles.edit_outersize,'String'));
FFTiv=handles.FFTiv;
FFTrv=handles.FFTrv;
IFFTiv=handles.IFFTiv;

cFrame = str2num(get(handles.edit_currentFrame,'String'));

switch handles.increment
    case 1
        endFrame   = str2num(get(handles.edit_lastframe,'String'));
    case -1
        endFrame   = str2num(get(handles.edit_firstframe,'String'));
end

if cFrame~=endFrame
    
    tp=cFrame;
    previousframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
    tp=cFrame+handles.increment;
    currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
    
    set(handles.edit_currentFrame,'String',num2str(tp));
    set(handles.slider_frames,'Value',tp);
    imshow(imadjust((currentframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
    
    opt = detbestlength2(FFTrv,FFTiv,IFFTiv,2*[outersize outersize],2*[cellsize cellsize],1,1);
    
    if  length(sisterList) >= tp && ~isempty(sisterList{tp}) && ~isempty(sisterList{1})
        sisExistInd = find(sisterList{tp}(:,1) ~= -1 & sisterList{tp}(:,1) ~= 0);
    else
        sisExistInd = [];
    end
    
    
    PosInd = find(cellpath{tp-handles.increment}(:,1)>0 & cellpath{tp-handles.increment}(:,2)>0)';
    
    for cell=PosInd
        if isempty(find(cell==sisExistInd,1))
            sisterList{tp}(cell,:) = sisterList{tp-handles.increment}(cell,:);
        end
        
        xL=max(cellpath{tp-handles.increment}(cell,1)-cellsize,1);
        xR=min(cellpath{tp-handles.increment}(cell,1)+cellsize,size(previousframe,2));
        yL=max(cellpath{tp-handles.increment}(cell,2)-cellsize,1);
        yR=min(cellpath{tp-handles.increment}(cell,2)+cellsize,size(previousframe,1));
        
        template = previousframe(yL:yR,xL:xR);
        
        xL=max(cellpath{tp-handles.increment}(cell,1)-outersize,1);
        xR=min(cellpath{tp-handles.increment}(cell,1)+outersize,size(currentframe,2));
        yL=max(cellpath{tp-handles.increment}(cell,2)-outersize,1);
        yR=min(cellpath{tp-handles.increment}(cell,2)+outersize,size(currentframe,1));
        
        testframe = currentframe(yL:yR,xL:xR);
        
        [x1 y1 maxVal] = corrMatching2(testframe, template,opt);
        if isempty(x1) || maxVal<str2num(get(handles.edit_similarityThres,'String'))
            x1 = cellpath{tp-handles.increment}(cell,1)-xL;
            y1 = cellpath{tp-handles.increment}(cell,2)-yL;
        end
        if get(handles.checkbox_autooptimize,'Value')==1
            [x1 y1 BW] = templateToCentroid(testframe,x1,y1,handles.maxI,handles.invertedLog);
        end
        
        cellpath{tp}(cell,:) = [xL+x1 yL+y1];
        assignin('base','cellpath',cellpath);
        
    end
    
    DeathInd = find(cellpath{tp-handles.increment}(:,1)==-2)';
    
    for cell=DeathInd
        sisterList{tp}(cell,:) = sisterList{tp-handles.increment}(cell,:);
        cellpath{tp}(cell,:) = cellpath{tp-handles.increment}(cell,:);
    end
    
    handles.sisterList = sisterList;
    handles.cellpath = cellpath;
    
    if ~isempty(bg)
        xshift = mean(cellpath{tp}(:,1))-mean(cellpath{tp-handles.increment}(:,1));
        yshift = mean(cellpath{tp}(:,2))-mean(cellpath{tp-handles.increment}(:,2));
        bg{tp}(:,1) = bg{tp-handles.increment}(:,1)+round(xshift);
        bg{tp}(:,2) = bg{tp-handles.increment}(:,2)+round(yshift);
        handles.bg = bg;
    end
    updateLists(cellpath,sisterList,res_cellpath,res_sisterList,bg,handles,tp);
    guidata(hObject, handles);
    handles = guidata(hObject);
    if get(handles.checkbox_cellmarking,'Value')
        [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
        handles.p = p;
        handles.bg_p = bg_p;
        guidata(hObject, handles);
    end
    drawnow;
    
    guidata(hObject, handles);
    handles = guidata(hObject);
    restoreCellList(handles,tp);
    
else
    set(handles.edit_commu,'String',{'Optimize points and press Record to save positions for each frame.'});
end



% --- Executes on button press in pushbutton_overlapframe.
function pushbutton_overlapframe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_overlapframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end
set(handles.togglebutton_editable,'Value',0);
cellpath = handles.cellpath;
sisterList = handles.sisterList;
res_cellpath = handles.res_cellpath;
res_sisterList = handles.res_sisterList;
bg = handles.bg;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));

cFrame = str2num(get(handles.edit_currentFrame,'String'));
switch handles.increment
    case 1
        endFrame   = str2num(get(handles.edit_lastframe,'String'));
    case -1
        endFrame   = str2num(get(handles.edit_firstframe,'String'));
end

if cFrame~=endFrame
    tp=cFrame+handles.increment;
    nextframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
    
    imshow(imadjust((nextframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
    set(handles.edit_currentFrame,'String',num2str(tp));
    set(handles.slider_frames,'Value',tp);
    
    if length(sisterList) >= tp && ~isempty(sisterList{tp}) && ~isempty(sisterList{1})
        sisExistInd = find(sisterList{tp}(:,1) ~= -1 & sisterList{tp}(:,1) ~= 0);
    else
        sisExistInd = [];
    end
    PosInd = find(cellpath{tp-handles.increment}(:,1)>0 & cellpath{tp-handles.increment}(:,2)>0)';
    
    for cell=PosInd
        if isempty(find(cell==sisExistInd,1))
            sisterList{tp}(cell,:) = sisterList{tp-handles.increment}(cell,:);
        end
        
        xL=max(cellpath{tp-handles.increment}(cell,1)-cellsize,1);
        xR=min(cellpath{tp-handles.increment}(cell,1)+cellsize,size(nextframe,2));
        yL=max(cellpath{tp-handles.increment}(cell,2)-cellsize,1);
        yR=min(cellpath{tp-handles.increment}(cell,2)+cellsize,size(nextframe,1));
        template = nextframe(yL:yR,xL:xR);
        switch handles.cellrecogmethod
            case 1
                [x1 y1 BW] = templateToCentroid(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog);
            case 2
                [x1 y1 BW] = templateToCentroid2(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog);
            case 3
                [x1 y1 BW] = templateToCentroid3(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog,str2num(get(handles.edit_outersize,'String')));
        end
        cellpath{tp}(cell,:) = [xL+x1 yL+y1];
        clear xL xR yL yR x1 x1 BW
    end
    
    DeathInd = find(cellpath{tp-handles.increment}(:,1)==-2)';
    
    for cell=DeathInd
        sisterList{tp}(cell,:) = sisterList{tp-handles.increment}(cell,:);
        cellpath{tp}(cell,:) = cellpath{tp-handles.increment}(cell,:);
    end
    
    handles.cellpath = cellpath;
    handles.sisterList = sisterList;
    
    if ~isempty(bg)
        xshift = mean(cellpath{tp}(:,1))-mean(cellpath{tp-handles.increment}(:,1));
        yshift = mean(cellpath{tp}(:,2))-mean(cellpath{tp-handles.increment}(:,2));
        bg{tp}(:,1) = bg{tp-handles.increment}(:,1)+round(xshift);
        bg{tp}(:,2) = bg{tp-handles.increment}(:,2)+round(yshift);
        handles.bg = bg;
    end
    guidata(hObject, handles);
    handles = guidata(hObject);
    if get(handles.checkbox_cellmarking,'Value')
        [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
        handles.p = p;
        handles.bg_p = bg_p;
        guidata(hObject, handles);
    end
    drawnow;
    updateLists(cellpath,sisterList,res_cellpath,res_sisterList,bg,handles,tp);
    guidata(hObject, handles);
    handles = guidata(hObject);
    restoreCellList(handles,tp);
else
    set(handles.edit_commu,'String',{'Optimize points and press Record to save positions for each frame.'});
end


% --- Executes on button press in pushbutton_overlapthru.
function pushbutton_overlapthru_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_overlapthru (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end
set(handles.togglebutton_editable,'Value',0);
cellpath = handles.cellpath;
sisterList = handles.sisterList;
res_cellpath = handles.res_cellpath;
res_sisterList = handles.res_sisterList;
bg = handles.bg;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));

cFrame = str2num(get(handles.edit_currentFrame,'String'));
switch handles.increment
    case 1
        endFrame   = str2num(get(handles.edit_lastframe,'String'))-1;
    case -1
        endFrame   = str2num(get(handles.edit_firstframe,'String'))+1;
end

if cFrame~=endFrame+handles.increment
    set(handles.edit_commu,'String',{'Processing...'});
    
    handles.greenflag = 1;
    guidata(hObject, handles);
    handles = guidata(hObject);
    
    for t=cFrame:handles.increment:endFrame
        tp=t+handles.increment;
        currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
        
        imshow(imadjust((currentframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
        set(handles.edit_currentFrame,'String',num2str(tp));
        set(handles.slider_frames,'Value',tp);
        
        
        if length(sisterList) >= tp && ~isempty(sisterList{tp}) && ~isempty(sisterList{1})
            sisExistInd = find(sisterList{tp}(:,1) ~= -1 & sisterList{tp}(:,1) ~= 0);
        else
            sisExistInd = [];
        end
        PosInd = find(cellpath{tp-handles.increment}(:,1)>0 & cellpath{tp-handles.increment}(:,2)>0)';
        
        for cell=PosInd
            if isempty(find(cell==sisExistInd,1))
                sisterList{tp}(cell,:) = sisterList{tp-handles.increment}(cell,:);
            end
            
            xL=max(cellpath{tp-handles.increment}(cell,1)-cellsize,1);
            xR=min(cellpath{tp-handles.increment}(cell,1)+cellsize,size(currentframe,2));
            yL=max(cellpath{tp-handles.increment}(cell,2)-cellsize,1);
            yR=min(cellpath{tp-handles.increment}(cell,2)+cellsize,size(currentframe,1));
            template = currentframe(yL:yR,xL:xR);
            switch handles.cellrecogmethod
                case 1
                    [x1 y1 BW] = templateToCentroid(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog);
                case 2
                    [x1 y1 BW] = templateToCentroid2(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog);
                case 3
                    [x1 y1 BW] = templateToCentroid3(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog,str2num(get(handles.edit_outersize,'String')));
            end
            cellpath{tp}(cell,:) = [xL+x1 yL+y1];
        end
        
        DeathInd = find(cellpath{tp-handles.increment}(:,1)==-2)';
        
        for cell=DeathInd
            sisterList{tp}(cell,:) = sisterList{tp-handles.increment}(cell,:);
            cellpath{tp}(cell,:) = cellpath{tp-handles.increment}(cell,:);
        end
        
        handles.cellpath = cellpath;
        handles.sisterList = sisterList;
        
        if ~isempty(bg)
            xshift = mean(cellpath{tp}(:,1))-mean(cellpath{tp-handles.increment}(:,1));
            yshift = mean(cellpath{tp}(:,2))-mean(cellpath{tp-handles.increment}(:,2));
            bg{tp}(:,1) = bg{tp-handles.increment}(:,1)+round(xshift);
            bg{tp}(:,2) = bg{tp-handles.increment}(:,2)+round(yshift);
            handles.bg = bg;
        end
        
        updateLists(cellpath,sisterList,res_cellpath,res_sisterList,bg,handles,tp);
        guidata(hObject, handles);
        handles = guidata(hObject);
        if get(handles.checkbox_cellmarking,'Value')
            [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
            handles.p = p;
            handles.bg_p = bg_p;
        end
        drawnow;
        
        handles = guidata(hObject);
        if handles.greenflag==0
            break;
        end
    end
    handles.greenflag = 1;
    
    
    if tp==endFrame
        set(handles.edit_commu,'String',{'Optimize points and press Record to save positions for each frame.'});
    end
    guidata(hObject, handles);
    handles = guidata(hObject);
    restoreCellList(handles,tp);
end
% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
first_tp   = str2num(get(handles.edit_firstframe,'String'));
last_tp   = str2num(get(handles.edit_lastframe,'String'));

H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
cellpath_name = ['/field' num2str(field) '/cellpath'];
sisterList_name = ['/field' num2str(field) '/sisterList'];
bg_name = ['/field' num2str(field) '/bg'];

SourceF = handles.SourceF;

if ~exist(fullfile(SourceF,H5filename),'file') 
    display([H5filename '.mat does not exist.']);
    return
else
    fileattrib(fullfile(SourceF,H5filename),'+w');
end

% Load Original seed points

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
    H5F.close(fid);
    cellpathinfo = h5info(fullfile(SourceF,H5filename), cellpath_name);
    
    cellpath_mat = h5read(fullfile(SourceF,H5filename),cellpath_name,[1 1 1], [cellpathinfo.Dataspace.Size(1) cellpathinfo.Dataspace.Size(2) cellpathinfo.Dataspace.Size(3)]);
    
    for tp=first_tp:cellpathinfo.Dataspace.Size(3)
        cellpath{tp} = cellpath_mat(:,:,tp);
    end
end

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,sisterList_name,'H5P_DEFAULT')
    H5F.close(fid);
    sisterListinfo = h5info(fullfile(SourceF,H5filename), sisterList_name);
    sisterList_mat = h5read(fullfile(SourceF,H5filename),sisterList_name,[1 1 1], [sisterListinfo.Dataspace.Size(1) sisterListinfo.Dataspace.Size(2) sisterListinfo.Dataspace.Size(3)]);
    
    for tp=first_tp:sisterListinfo.Dataspace.Size(3)
        sisterList{tp} = sisterList_mat(:,:,tp);
    end
end

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,bg_name,'H5P_DEFAULT')
    H5F.close(fid);
    bginfo = h5info(fullfile(SourceF,H5filename), bg_name);
    bg_mat = h5read(fullfile(SourceF,H5filename),bg_name,[1 1 1], [bginfo.Dataspace.Size(1) bginfo.Dataspace.Size(2) bginfo.Dataspace.Size(3)]);
    
    for tp=first_tp:bginfo.Dataspace.Size(3)
        bg{tp} = bg_mat(:,:,tp);
    end
end


handles.cellpath=cellpath;
handles.sisterList=sisterList;
handles.bg = bg;

tp = str2num(get(handles.edit_currentFrame,'String'));
set(handles.togglebutton_editable,'Value',0);

[p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,1);
handles.p = p;
handles.bg_p = bg_p;
set(handles.edit_commu,'String',['Loaded data from ' H5filename]);
guidata(hObject, handles);
handles = guidata(hObject);
set(handles.listbox_cells,'Value',1);
updateLists(cellpath,sisterList,[],[],bg,handles,tp)
restoreCellList(handles,tp);


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg = handles.bg;

first_tp   = str2num(get(handles.edit_firstframe,'String'));
last_tp   = str2num(get(handles.edit_lastframe,'String'));
c_tp = str2num(get(handles.edit_currentFrame,'String'));

cellpath_mat = -1*(ones(size(cellpath{c_tp},1),2,length(cellpath)));
sisterList_mat = -1*(ones(size(sisterList{c_tp},1),size(sisterList{c_tp},2),length(sisterList)));
bg_mat = -1*(ones(size(bg{c_tp},1),2,length(bg)));

for tp=1:length(cellpath)
    if ~isempty(cellpath{tp})
        cellpath_mat(:,:,tp) = cellpath{tp};
        sisterList_mat(:,:,tp) = sisterList{tp};
        bg_mat(:,:,tp) = bg{tp};
    end
end

H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
cellpath_name = ['/field' num2str(field) '/cellpath'];
sisterList_name = ['/field' num2str(field) '/sisterList'];
bg_name = ['/field' num2str(field) '/bg'];

if exist(fullfile(handles.SourceF,H5filename),'file')
    fileattrib(fullfile(handles.SourceF,H5filename),'+w');
    
    fid = H5F.open(fullfile(handles.SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if ~H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
        H5F.close(fid);
        display(['Initializing ' H5filename ':' cellpath_name]);
    else
        H5L.delete(fid,cellpath_name,'H5P_DEFAULT');
        display(['Overwriting ' H5filename ':' cellpath_name]);
        H5F.close(fid);
    end
end

h5create(fullfile(handles.SourceF,H5filename), cellpath_name, [size(cellpath_mat,1), size(cellpath_mat,2), size(cellpath_mat,3)], 'Datatype', 'double', 'ChunkSize', [1, size(cellpath_mat,2), size(cellpath_mat,3)], 'Deflate', 9);
h5write(fullfile(handles.SourceF,H5filename), cellpath_name, cellpath_mat, [1 1 1], [size(cellpath_mat,1) size(cellpath_mat,2) size(cellpath_mat,3)]);

fid = H5F.open(fullfile(handles.SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,sisterList_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' sisterList_name]);
else
    H5L.delete(fid,sisterList_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' sisterList_name]);
    H5F.close(fid);
end

h5create(fullfile(handles.SourceF,H5filename), sisterList_name, [size(sisterList_mat,1), size(sisterList_mat,2), size(sisterList_mat,3)], 'Datatype', 'double', 'ChunkSize', [1, size(sisterList_mat,2), size(sisterList_mat,3)], 'Deflate', 9);
h5write(fullfile(handles.SourceF,H5filename), sisterList_name, sisterList_mat, [1 1 1], [size(sisterList_mat,1) size(sisterList_mat,2) size(sisterList_mat,3)]);

fid = H5F.open(fullfile(handles.SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,bg_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' bg_name]);
else
    H5L.delete(fid,bg_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' bg_name]);
    H5F.close(fid);
end

h5create(fullfile(handles.SourceF,H5filename), bg_name, [size(bg_mat,1), size(bg_mat,2), size(bg_mat,3)], 'Datatype', 'double', 'ChunkSize', [1, size(bg_mat,2), size(bg_mat,3)], 'Deflate', 9);
h5write(fullfile(handles.SourceF,H5filename), bg_name, bg_mat, [1 1 1], [size(bg_mat,1) size(bg_mat,2) size(bg_mat,3)]);

set(handles.edit_commu,'String',['Your data is saved to ' H5filename]);

function edit_cellNo_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cellNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cellNo as text
%        str2double(get(hObject,'String')) returns contents of edit_cellNo as a double
p=handles.p;
cellpath = handles.cellpath;
sisterList = handles.sisterList;
res_cellpath = handles.res_cellpath;
res_sisterList = handles.res_sisterList;

tp = str2num(get(handles.edit_currentFrame,'String'));
cell = str2num(get(hObject,'String'));

if isempty(cellpath) || isempty(cellpath{tp})
    set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',res_cellpath{tp}(cell,1),res_cellpath{tp}(cell,2))));
    set(handles.edit_sister,'String',num2str(res_sisterList{tp}(cell,:)));
    set(handles.listbox_Restcells,'Value',cell);
    
else
    set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',cellpath{tp}(cell,1),cellpath{tp}(cell,2))));
    set(handles.edit_sister,'String',num2str(sisterList{tp}(cell,:)));
    set(handles.listbox_cells,'Value',cell);
end

guidata(hObject, handles);

% --- Executes on button press in pushbutton_gensister.
function pushbutton_gensister_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_gensister (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.togglebutton_editable,'Value',1);
[p bg_p] = changeToEditable(handles);
handles.p = p;
handles.bg_p = bg_p;
guidata(hObject, handles);
handles = guidata(hObject);
p = handles.p;
cellpath = handles.cellpath;
sisterList = handles.sisterList;
c_tp = str2num(get(handles.edit_currentFrame,'String'));
cell   = str2num(get(handles.edit_cellNo,'String'));
sister = str2num(get(handles.edit_sister,'String'));

first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));
switch handles.increment
    case 1
        for tp=first_tp:1:length(cellpath)
            if ~isempty(cellpath{tp})
                endFrame = tp;
            end
        end
    case -1
        endFrame   = str2num(get(handles.edit_firstframe,'String'));
        for tp=last_tp:-1:first_tp
            if ~isempty(cellpath{tp})
                endFrame = tp;
            end
        end
end


totalSize = size(cellpath{tp},1);
currentSize = size(cellpath{c_tp},1);
[xi,yi,but] = ginput(1);
%add new points to the list
newcell = totalSize+1;
if newcell > currentSize+1
    for tp = first_tp:endFrame
        for c = (currentSize+1):totalSize
            if size(cellpath{tp},1) < c
                cellpath{tp}(c,:) =[-1 -1];
                sisterList{tp}(c,:)    = [-1 -1 -1];
            end
        end
    end
end



cellpath{c_tp}(newcell,:) = round([xi yi]);
for tp = first_tp:c_tp-1
    cellpath{tp}(newcell,:) =[-1 -1];
    sisterList{tp}(newcell,:)    = [-1 -1 -1];
end

%check if 'cell' already has a sister


if sisterList{c_tp}(cell,1)==-1
    %update sisterList by pointing to each other
    sisterList{c_tp}(cell,1) = newcell;
    sisterList{c_tp}(newcell,1) = cell;
    sisterList{c_tp}(newcell,2) = -1;
    sisterList{c_tp}(newcell,3) = -1;
elseif sisterList{c_tp}(cell,2)==-1
    %update sisterList by pointing to each other
    sisterList{c_tp}(cell,2)    = newcell;
    sisterList{c_tp}(newcell,1) = sisterList{c_tp}(cell,1);
    sisterList{c_tp}(newcell,2) = cell;
    sisterList{c_tp}(newcell,3) = -1;
elseif sisterList{c_tp}(cell,3)==-1
    %update sisterList by pointing to each other
    sisterList{c_tp}(cell,3)    = newcell;
    sisterList{c_tp}(newcell,1) = sisterList{c_tp}(cell,1);
    sisterList{c_tp}(newcell,2) = sisterList{c_tp}(cell,2);
    sisterList{c_tp}(newcell,3) = cell;
else
    set(handles.edit_commu,'String','ERROR:No space to store 4th-gen sister.');
end

set(handles.edit_sister,'String',num2str(sisterList{c_tp}(cell,:)));
p{newcell} = coordtopoint(cellpath{c_tp}(newcell,:),newcell,'r',handles.fcn1,handles.axes1,get(handles.checkbox_cellnostring,'Value'));
addNewPositionCallback(p{newcell},@(h) set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',h(1),h(2)))));
addNewPositionCallback(p{newcell},@(h) set(handles.edit_cellNo,'String',num2str(newcell)));
addNewPositionCallback(p{newcell},@(h) set(handles.listbox_cells,'Value',newcell));
addNewPositionCallback(p{newcell},@(h) set(handles.edit_sister,'String',num2str(sisterList{c_tp}(newcell,:))));
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,handles.bg,handles,c_tp);
handles.p = p;
handles.cellpath = cellpath;
handles.sisterList = sisterList;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listbox_cells_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to listbox_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
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



function edit_CH_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CH as text
%        str2double(get(hObject,'String')) returns contents of edit_CH as a double


% --- Executes during object creation, after setting all properties.
function edit_CH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CH (see GCBO)
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


function edit_firstframe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_firstframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_firstframe as text
%        str2double(get(hObject,'String')) returns contents of edit_firstframe as a double

minF = str2num(get(hObject,'String'));
maxF = str2num(get(handles.edit_lastframe,'String'));
set(handles.slider_frames,'Min',minF);
set(handles.slider_frames,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);

% --- Executes during object creation, after setting all properties.
function edit_firstframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_firstframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_lastframe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lastframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lastframe as text
%        str2double(get(hObject,'String')) returns contents of
%        edit_lastframe as a double.
maxF = str2num(get(hObject,'String'));
minF = str2num(get(handles.edit_firstframe,'String'));
set(handles.slider_frames,'Max',maxF);
set(handles.slider_frames,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);




% --- Executes during object creation, after setting all properties.
function edit_lastframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lastframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_thresLow_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thresLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thresLow as text
%        str2double(get(hObject,'String')) returns contents of edit_thresLow as a double

handles.thresLow = str2num(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_thresLow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thresLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.thresLow = str2num(get(hObject,'String'));
guidata(hObject, handles);



function edit_thresHigh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thresHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thresHigh as text
%        str2double(get(hObject,'String')) returns contents of edit_thresHigh as a double

handles.thresHigh = str2num(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_thresHigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thresHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.thresHigh = str2num(get(hObject,'String'));
guidata(hObject, handles);



function edit_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sigma as text
%        str2double(get(hObject,'String')) returns contents of edit_sigma as a double
handles.sigma = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.sigma = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes on button press in checkbox_inverting.
function checkbox_inverting_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_inverting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_inverting
handles.invertedLog = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function checkbox_inverting_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_inverting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.invertedLog = get(hObject,'Value');
guidata(hObject, handles);


function edit_thresMax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thresMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thresMax as text
%        str2double(get(hObject,'String')) returns contents of edit_thresMax as a double

cellpath = handles.cellpath;
bg = handles.bg;
sisterList = handles.sisterList;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
c_tp = str2num(get(handles.edit_currentFrame,'String'));
set(handles.slider_frames,'Value',str2num(get(handles.edit_currentFrame,'String')));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));
if c_tp>last_tp
    tp=last_tp;
elseif c_tp<first_tp
    tp=first_tp;
else
    tp=c_tp;
end

previousframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
set(handles.edit_currentFrame,'String',num2str(tp));
imshow(imadjust((previousframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(hObject,'String'))],[0 1]),'Parent',handles.axes1);
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;

% --- Executes during object creation, after setting all properties.
function edit_thresMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thresMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_cellareathres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cellareathres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.cellareathres = str2num(get(hObject,'String'));
guidata(hObject, handles);



function edit_thresstep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thresstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thresstep as text
%        str2double(get(hObject,'String')) returns contents of edit_thresstep as a double
handles.thresstep = str2num(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_thresstep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thresstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.thresstep = str2num(get(hObject,'String'));
guidata(hObject, handles);

function edit_cellareathres_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cellareathres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cellareathres as text
%        str2double(get(hObject,'String')) returns contents of edit_cellareathres as a double
handles.cellareathres = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_currentFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_currentFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function edit_cellNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cellNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_coord_Callback(hObject, eventdata, handles)
% hObject    handle to edit_coord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_coord as text
%        str2double(get(hObject,'String')) returns contents of edit_coord as a double
cellpath = handles.cellpath;
tp = str2num(get(handles.edit_currentFrame,'String'));
selected_cell = str2num(get(handles.edit_cellNo,'String'));

currentCoord = get(hObject,'String');
tokens   = regexp(currentCoord, '((?<myX>\d+),(?<myY>\d+))|(?<myX>\d+),(?<myY>\d+)','tokens');
myX = str2num(tokens{1}{1});
myY = str2num(tokens{1}{2});
cellpath{tp}(selected_cell,:) = [myX myY];
set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',myX,myY)));
handles.cellpath = cellpath;
guidata(hObject, handles);
handles = guidata(hObject);
restoreCellList(handles,tp);
updateLists(cellpath,handles.sisterList,handles.res_cellpath,handles.res_sisterList,handles.bg,handles,tp);

% --- Executes during object creation, after setting all properties.
function edit_coord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_coord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_FRETCH_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FRETCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FRETCH as text
%        str2double(get(hObject,'String')) returns contents of edit_FRETCH as a double


% --- Executes during object creation, after setting all properties.
function edit_FRETCH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FRETCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_CFPCH_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CFPCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CFPCH as text
%        str2double(get(hObject,'String')) returns contents of edit_CFPCH as a double


% --- Executes during object creation, after setting all properties.
function edit_CFPCH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CFPCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sister_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sister (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sister as text
%        str2double(get(hObject,'String')) returns contents of edit_sister as a double

p=handles.p;
cellpath = handles.cellpath;
sisterList = handles.sisterList;
tp = str2num(get(handles.edit_currentFrame,'String'));
selected_cell = str2num(get(handles.edit_cellNo,'String'));
sisterList{tp}(selected_cell,:) = str2num(get(hObject,'String'));
handles.sisterList = sisterList;
guidata(hObject, handles);
handles = guidata(hObject);
restoreCellList(handles,tp);
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,handles.bg,handles,tp);
% --- Executes during object creation, after setting all properties.
function edit_sister_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sister (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cellsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cellsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cellsize as text
%        str2double(get(hObject,'String')) returns contents of edit_cellsize as a double


% --- Executes during object creation, after setting all properties.
function edit_cellsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cellsize (see GCBO)
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



function edit_fileformat_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fileformat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fileformat as text
%        str2double(get(hObject,'String')) returns contents of edit_fileformat as a double


% --- Executes during object creation, after setting all properties.
function edit_fileformat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fileformat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_cellnostring.
function checkbox_cellnostring_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cellnostring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_cellnostring


% --- Executes on button press in pushbutton_genim.
function pushbutton_genim_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_genim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end

set(handles.edit_commu,'String',{'Generating printable plot of current frame...'});
cellpath = handles.cellpath;
sisterList = handles.sisterList;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
tp = str2num(get(handles.edit_currentFrame,'String'));
fileformat = get(handles.edit_fileformat,'String');

first_tp = str2num(get(handles.edit_firstframe,'String'));
if handles.filetype == 3
    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},first_tp);
    else
        filename = sprintf(fileformat,first_tp);
    end
    first_info = imfinfo(fullfile(handles.SourceF,filename));
end

currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
if handles.filetype == 3
    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},tp);
    else
        filename = sprintf(fileformat,tp);
    end
    current_info = imfinfo(fullfile(handles.SourceF,filename));
end

figure(1),imshow(imadjust((currentframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]));hold on;
rancolor = lines(24);
if ~isempty(cellpath)
    for cell=1:size(cellpath{tp},1)
        if sisterList{tp}(cell)==-1
            plot(cellpath{tp}(cell,1),cellpath{tp}(cell,2),'rx');
        else
            plot(cellpath{tp}(cell,1),cellpath{tp}(cell,2),'o','MarkerFaceColor',rancolor(mod(cell,24),:),'MarkerEdgeColor',rancolor(mod(cell,24),:));
            plot(cellpath{tp}(sisterList{tp}(cell),1),cellpath{tp}(sisterList{tp}(cell),2),'o','MarkerFaceColor',rancolor(mod(cell,24),:),'MarkerEdgeColor',rancolor(mod(cell,24),:));
        end
        text(cellpath{tp}(cell,1)+10,cellpath{tp}(cell,2)+10,num2str(cell),'HorizontalAlignment','left',...
            'VerticalAlignment','middle','color',[0 .9 .5]);
    end
end
switch handles.framestamp
    case 1
        text(30,30,['frame:' num2str(tp)],...
            'HorizontalAlignment','left','VerticalAlignment','top','color','w','BackgroundColor','k','fontsize',22);
    case 2
        switch handles.filetype
            case 3
                [Y, M, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
                hour = 24*D+round(H);
                minute = round(MN);
                second = round(S);
            otherwise
                hour = floor(1/60*(tp-1)*( str2num(get(handles.edit_minstamp,'String')) + str2num(get(handles.edit_secstamp,'String'))/60 ));
                minute = mod(floor((tp-1)*( str2num(get(handles.edit_minstamp,'String')) + str2num(get(handles.edit_secstamp,'String'))/60 )),60);
                second = mod((tp-1)*( str2num(get(handles.edit_minstamp,'String'))*60 + str2num(get(handles.edit_secstamp,'String'))),60);
        end
        text(30,30,[num2str(hour,'%02.0f') ':' num2str(minute,'%02.0f') ':' num2str(second,'%02.0f')],...
            'HorizontalAlignment','left','VerticalAlignment','top','color','w','BackgroundColor','k','fontsize',22);
end

set(handles.edit_commu,'String',{'Use menu bar of the new figure to save plot'});

% --- Executes on button press in pushbutton_genmov.
function pushbutton_genmov_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_genmov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end

set(handles.edit_commu,'String',{'Generating MATLAB movie...'});
cellpath = handles.cellpath;
sisterList = handles.sisterList;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
fileformat = get(handles.edit_fileformat,'String');
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));
n=1;

aviobj = VideoWriter(fullfile(get(handles.edit_sourceF,'String'),['myMov_r' num2str(row) 'c' num2str(col) 'f' num2str(field) 'ch' num2str(channel) '.avi']));
aviobj.FrameRate = 12;
aviobj.Quality = 80;
open(aviobj);
f1=figure();
axes(gca);
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

if handles.filetype == 3
    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},first_tp);
    else
        filename = sprintf(fileformat,first_tp);
    end
    first_info = imfinfo(fullfile(handles.SourceF,filename));
end

for tp=first_tp:last_tp
    currentframe = loadimage(handles.filetype,fileformat,[row col field plane channel],tp,handles.channelnames,handles.SourceF);
    if handles.filetype == 3
        if ~isempty(handles.channelnames)
            filename = sprintf(fileformat,handles.channelnames{channel},tp);
        else
            filename = sprintf(fileformat,tp);
        end
        current_info = imfinfo(fullfile(handles.SourceF,filename));
    end
    
    figure(1),imshow(imadjust((currentframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'InitialMagnification',50);hold on;
    if ~isempty(cellpath)
        for cell=1:size(cellpath{tp},1)
            if sisterList{tp}(cell)==-1
                plot(cellpath{tp}(cell,1),cellpath{tp}(cell,2),'rx');
            else
                plot(cellpath{tp}(cell,1),cellpath{tp}(cell,2),'go');
            end
            text(cellpath{tp}(cell,1)+10,cellpath{tp}(cell,2)+10,num2str(cell),'HorizontalAlignment','left',...
                'VerticalAlignment','middle','color',[0 .9 .5]);
        end
    end
    switch handles.framestamp
        case 1
            text(30,30,['frame:' num2str(tp)],...
                'HorizontalAlignment','left','VerticalAlignment','top','color','w','BackgroundColor','k','fontsize',22);
        case 2
            switch handles.filetype
                case 3
                    [Y, M, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
                    hour = 24*D+round(H);
                    minute = round(MN);
                    second = round(S);
                otherwise
                    hour = floor(1/60*(tp-1)*( str2num(get(handles.edit_minstamp,'String')) + str2num(get(handles.edit_secstamp,'String'))/60 ));
                    minute = mod(floor((tp-1)*( str2num(get(handles.edit_minstamp,'String')) + str2num(get(handles.edit_secstamp,'String'))/60 )),60);
                    second = mod((tp-1)*( str2num(get(handles.edit_minstamp,'String'))*60 + str2num(get(handles.edit_secstamp,'String'))),60);
            end
            text(30,30,[num2str(hour,'%02.0f') ':' num2str(minute,'%02.0f') ':' num2str(second,'%02.0f')],...
                'HorizontalAlignment','left','VerticalAlignment','top','color','w','BackgroundColor','k','fontsize',22);
    end
    drawnow;
    F = getframe(gca);
    writeVideo(aviobj,F);
end

close(aviobj);
set(handles.edit_commu,'String',['Done generating: myMov_r' num2str(row) 'c' num2str(col) 'f' num2str(col) 'ch' num2str(channel) '.avi']);
close(f1);
% --- Executes during object creation, after setting all properties.
function edit_movecoord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_movecoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_outersize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_outersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxcells_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxcells as text
%        str2double(get(hObject,'String')) returns contents of edit_maxcells as a double


% --- Executes during object creation, after setting all properties.
function edit_maxcells_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_randomized.
function checkbox_randomized_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_randomized (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_randomized


% --- Executes when selected object is changed in uibuttongroup_filetype.
function uibuttongroup_filetype_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup_filetype
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton_petiffs'
        set(handles.edit_commu,'String','PE tiffs chosen for inputs');
        set(handles.edit_fileformat,'String','r%02.0fc%02.0ff%02.0fp%02.0frc%1.0f-ch1sk%ufk1fl1.tiff');
        handles.filetype = 1;
    case 'radiobutton_tiffstack'
        set(handles.edit_commu,'String','DV stacked tiff chosen for input. Please locate your file.');
        handles.filetype = 2;
    case 'radiobutton_customtiff'
        set(handles.edit_commu,'String','Metamorph tiff chosen for inputs. Edit file structure if needed.');
        handles.filetype = 3;
    case 'radiobutton_columbus'
        set(handles.edit_commu,'String','Columbus tif chosen for inputs. Edit file structure if needed.');
        set(handles.edit_fileformat,'String','%03.0f%03.0f-%u-%03.0f%03.0f%03.0f.tif');
        handles.filetype = 4;     
        
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton_locatefile.
function pushbutton_locatefile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_locatefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname,FilterIndex] = uigetfile(...
    {'*.tif;*.fig;*.mat;*.mdl','stacked tiff Files (*.tif)';
    '*.*',  'All Files (*.*)'}, 'Pick a file');
if FilterIndex~=0
    set(handles.edit_fileformat,'String',filename);
    handles.tiffstackname = fullfile(pathname, filename);
end
guidata(hObject, handles);



function edit_totalCHs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_totalCHs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_totalCHs as text
%        str2double(get(hObject,'String')) returns contents of edit_totalCHs as a double


% --- Executes during object creation, after setting all properties.
function edit_totalCHs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_totalCHs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uibuttongroup_filetype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup_filetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function uibuttongroup_cellrecogmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup_cellrecogmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.cellrecogmethod =1 ;
guidata(hObject, handles);



function edit_outersize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_outersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_outersize as text
%        str2double(get(hObject,'String')) returns contents of edit_outersize as a double


% --- Executes on selection change in listbox_bgs.
function listbox_bgs_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_bgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_bgs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_bgs

% --- Executes during object creation, after setting all properties.
function listbox_bgs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_bgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton_bgadd.
function pushbutton_bgadd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_bgadd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bg = handles.bg;
bg_p = handles.bg_p;
tp = str2num(get(handles.edit_currentFrame,'String'));
if isempty(bg)
    n=0;
else
    n=size(bg{tp},1);
end

[xi,yi] = ginput(1);
bg{tp}(n+1,:) = round([xi yi]);
bg_p{n+1} = coordtopoint([xi yi],n+1,'g',handles.fcn1,handles.axes1,get(handles.checkbox_cellnostring,'Value'));
set(handles.listbox_bgs,'String',num2str( [(1:size(bg{tp},1))' bg{tp}] ),'Value',n+1);
addNewPositionCallback(bg_p{n+1},@(h) set(handles.listbox_bgs,'Value',n+1));
handles.bg = bg;
handles.bg_p = bg_p;
guidata(hObject, handles);

% --- Executes on button press in pushbutton_bgremove.
function pushbutton_bgremove_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_bgremove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bg_p=handles.bg_p;
bg = handles.bg;
tp = str2num(get(handles.edit_firstframe,'String'));
del_b = get(handles.listbox_bgs,'Value');

oldsize = length(bg{tp});
if del_b<oldsize
    temp(1:(del_b-1),:) = bg{tp}(1:(del_b-1),:);
    temp(del_b:(oldsize-1),:) = bg{tp}((del_b+1):oldsize,:);
else
    temp(1:(del_b-1),:) = bg{tp}(1:(del_b-1),:);
    
end
bg{tp} = temp;

for b=1:length(bg_p)
    delete(bg_p{b});
end
clear bg_p;

if ~isempty(bg)
    set(handles.listbox_bgs,'String',num2str( [(1:size(bg{tp},1))' bg{tp}] ));
    for b=1:size(bg{tp},1)
        bg_p{b} = coordtopoint(bg{tp}(b,:),b,'g',handles.fcn1,handles.axes1,get(handles.checkbox_cellnostring,'Value'));
    end
    
end
set(handles.listbox_bgs,'Value',1);
handles.bg = bg;
handles.bg_p = bg_p;
guidata(hObject, handles);


% --- Executes on slider movement.
function slider_frames_Callback(hObject, eventdata, handles)
% hObject    handle to slider_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.edit_currentFrame,'String',num2str(round(get(hObject,'Value'))));

cellpath = handles.cellpath;
bg = handles.bg;
sisterList = handles.sisterList;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
c_tp = str2num(get(handles.edit_currentFrame,'String'));
fileformat = get(handles.edit_fileformat,'String');
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

if handles.filetype == 3
    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},first_tp);
    else
        filename = sprintf(fileformat,first_tp);
    end
    first_info = imfinfo(fullfile(handles.SourceF,filename));
end
if c_tp>last_tp
    tp=last_tp;
elseif c_tp<first_tp
    tp=first_tp;
else
    tp=c_tp;
end

currentframe = loadimage(handles.filetype,fileformat,[row col field plane channel],tp,handles.channelnames,handles.SourceF);
if handles.filetype == 3
    if ~isempty(handles.channelnames)
        filename = sprintf(fileformat,handles.channelnames{channel},tp);
    else
        filename = sprintf(fileformat,tp);
    end
    current_info = imfinfo(fullfile(handles.SourceF,filename));
    [Y, M, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
    hour = 24*D+round(H);
    minute = round(MN);
    second = round(S);
    set(handles.edit_commu,'String',[num2str(hour,'%02.0f') ':' num2str(minute,'%02.0f') ':' num2str(second,'%02.0f')]);
end
imshow(imadjust((currentframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
set(handles.edit_currentFrame,'String',num2str(tp));

if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,tp);

% --- Executes during object creation, after setting all properties.
function slider_frames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes during object creation, after setting all properties.
function CellTracking_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellTracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox_autooptimize.
function checkbox_autooptimize_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_autooptimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_autooptimize



function edit_similarityThres_Callback(hObject, eventdata, handles)
% hObject    handle to edit_similarityThres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_similarityThres as text
%        str2double(get(hObject,'String')) returns contents of edit_similarityThres as a double


% --- Executes during object creation, after setting all properties.
function edit_similarityThres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_similarityThres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_cellmarking.
function checkbox_cellmarking_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cellmarking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_cellmarking


% --- Executes on button press in togglebutton_backward.
function togglebutton_backward_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_backward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg = handles.bg;

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));

currentFrame = str2num(get(handles.edit_currentFrame,'String'));
firstframe   = str2num(get(handles.edit_firstframe,'String'));

nextframe = currentFrame-1;
if nextframe<firstframe
    nextframe=firstframe;
end
if nextframe == currentFrame
    return;
end

FFTiv=handles.FFTiv;
FFTrv=handles.FFTrv;
IFFTiv=handles.IFFTiv;

previousframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],currentFrame,handles.channelnames,handles.SourceF);
currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],nextframe,handles.channelnames,handles.SourceF);

opt = detbestlength2(FFTrv,FFTiv,IFFTiv,size(imresize(currentframe,0.1)),size(imresize(previousframe,0.1)),isreal(imresize(currentframe,0.1)),isreal(imresize(previousframe,0.1)));
[x y maxVal] = corrMatching2(imresize(currentframe,0.1),imresize(previousframe,0.1),opt);

xshift = 10*x-round(size(currentframe,2)/2);
yshift = 10*y-round(size(currentframe,1)/2);

cellpath{nextframe} = [cellpath{currentFrame}(:,1)+xshift cellpath{currentFrame}(:,2)+yshift];
sisterList{nextframe} = sisterList{currentFrame};

if ~isempty(bg)
    bg{nextframe} = [bg{currentFrame}(:,1)+xshift bg{currentFrame}(:,2)+yshift];
    set(handles.listbox_bgs,'String',num2str( [(1:size(bg{nextframe},1))' bg{nextframe}] ));
end
set(handles.listbox_cells,'String',num2str( [(1:size(cellpath{nextframe},1))' cellpath{nextframe} sisterList{nextframe}] ));
handles.cellpath = cellpath;
handles.sisterList = sisterList;
handles.bg = bg;
guidata(hObject, handles);
handles = guidata(hObject);
restoreCellList(handles,nextframe);
imshow(imadjust((currentframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
set(handles.edit_currentFrame,'String',num2str(nextframe));
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,nextframe,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;

% --- Executes on button press in togglebutton_editable.
function togglebutton_editable_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_editable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end
[p bg_p] = changeToEditable(handles);
handles.p = p;
handles.bg_p = bg_p;
guidata(hObject, handles);



function [p bg_p] = changeToEditable(handles)
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
tp = str2num(get(handles.edit_currentFrame,'String'));
cellpath = handles.cellpath;
bg = handles.bg;
sisterList = handles.sisterList;
currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
imshow(imadjust((currentframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
    drawnow;
end



% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_thresMin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thresMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thresMin as text
%        str2double(get(hObject,'String')) returns contents of edit_thresMin as a double
cellpath = handles.cellpath;
bg = handles.bg;
sisterList = handles.sisterList;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
c_tp = str2num(get(handles.edit_currentFrame,'String'));
set(handles.slider_frames,'Value',str2num(get(handles.edit_currentFrame,'String')));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));
if c_tp>last_tp
    tp=last_tp;
elseif c_tp<first_tp
    tp=first_tp;
else
    tp=c_tp;
end

previousframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
set(handles.edit_currentFrame,'String',num2str(tp));
imshow(imadjust((previousframe),[str2num(get(hObject,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;

% --- Executes during object creation, after setting all properties.
function edit_thresMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thresMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

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
handles.maxI = str2num(get(hObject,'String'));
guidata(hObject, handles);

function edit_maxI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxI as text
%        str2double(get(hObject,'String')) returns contents of edit_maxI as a double
handles.maxI = str2num(get(hObject,'String'));
guidata(hObject, handles);



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
        
        testInd = regexp(tline,'WaveName\d+');
        num = length(testInd);
        if num > 0
            wavename1  = regexp(tline, '(?<="WaveName\d+", ")\w+(?=_)', 'match');
            wavename2  = regexp(tline, '(?<="WaveName\d+", "\w+_).+(?=")', 'match');
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

[notp stagePos stageName waveName] = readndfile(fullfile(handles.SourceF,handles.ndfilename));

if notp==-1
    return;
end


handles.filetype = 3;
set(handles.radiobutton_customtiff,'Value',1);
set(handles.edit_firstframe,'String',num2str(1));
set(handles.edit_lastframe,'String',num2str(notp));
set(handles.edit_currentFrame,'String',num2str(1));
set(handles.popupmenu_stagePos,'String',stagePos);
set(handles.popupmenu_stagePos,'Value',1);
set(handles.edit_stageInfo,'String',stageName{1});
handles.stageName = stageName;
handles.channelnames = waveName;

prefix = handles.ndfilename(1:(end-3));
handles.prefix = prefix;

if ~isempty(waveName)
    fileformat = [prefix '_%s_s1_t%g.TIF'];
else
    fileformat = [prefix '_s1_t%g.TIF'];
end
set(handles.edit_fileformat,'String',fileformat);

L = regexp(stageName{1}, 'r(?<row>\d+)','names');
if ~isempty(L)
    row = L.row;
    set(handles.edit_row,'String',row);
else
    set(handles.edit_row,'String','1');
end

L = regexp(stageName{1}, 'c(?<col>\d+)','names');
if ~isempty(L)
    col = L.col;
    set(handles.edit_col,'String',col);
else
    set(handles.edit_col,'String','1');
end
L = regexp(stageName{1}, 'f(?<field>\d+)','names');

if ~isempty(L)
    field = L.field;
    set(handles.edit_field,'String',field);
else
    set(handles.edit_field,'String','1');
end

guidata(hObject, handles);


% --- Executes on selection change in popupmenu_stagePos.
function popupmenu_stagePos_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_stagePos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.edit_stageInfo,'String',handles.stageName{get(hObject,'Value')});


if ~isempty(handles.channelnames)
    fileformat = [handles.prefix '_%s_s' num2str(get(hObject,'Value')) '_t%g.TIF'];
else
    fileformat = [handles.prefix '_s' num2str(get(hObject,'Value')) '_t%g.TIF'];
end

set(handles.edit_fileformat,'String',fileformat);

L = regexp(handles.stageName{get(hObject,'Value')}, 'r(?<row>\d+)','names');

if ~isempty(L)
    row = L.row;
    set(handles.edit_row,'String',row);
else
    set(handles.edit_row,'String','1');
end

L = regexp(handles.stageName{get(hObject,'Value')}, 'c(?<col>\d+)','names');

if ~isempty(L)
    col = L.col;
    set(handles.edit_col,'String',col);
else
    set(handles.edit_col,'String','1');
end
L = regexp(handles.stageName{get(hObject,'Value')}, 'f(?<field>\d+)','names');

if ~isempty(L)
    field = L.field;
    set(handles.edit_field,'String',field);
else
    set(handles.edit_field,'String','1');
end





guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_stagePos contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_stagePos


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


% --- Executes on button press in pushbutton_duplicateFbF.
function pushbutton_duplicateFbF_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_duplicateFbF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end
set(handles.togglebutton_editable,'Value',0);
cellpath = handles.cellpath;
sisterList = handles.sisterList;
res_cellpath = handles.res_cellpath;
res_sisterList = handles.res_sisterList;

bg = handles.bg;

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));

currentFrame = str2num(get(handles.edit_currentFrame,'String'));
firstframe = str2num(get(handles.edit_firstframe,'String'));
lastframe   = str2num(get(handles.edit_lastframe,'String'));

nextframe = currentFrame+handles.increment;

switch handles.increment
    case 1
        if nextframe>lastframe
            nextframe=lastframe;
        end
    case -1
        if nextframe<firstframe
            nextframe=firstframe;
        end
end
if nextframe == currentFrame
    return;
end

switch handles.increment
    case 1
        endframe=lastframe-1;
    case -1
        endframe=firstframe+1;
end




FFTiv=handles.FFTiv;
FFTrv=handles.FFTrv;
IFFTiv=handles.IFFTiv;


previousframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],currentFrame,handles.channelnames,handles.SourceF);
previousframe = previousframe(60:end-60,60:end-60);
currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],nextframe,handles.channelnames,handles.SourceF);

opt = detbestlength2(FFTrv,FFTiv,IFFTiv,size(imresize(currentframe,1)),size(imresize(previousframe,1)),isreal(imresize(currentframe,1)),isreal(imresize(previousframe,1)));
[x y maxVal] = corrMatching2(imresize(currentframe,1),imresize(previousframe,1),opt);

xshift = x-round(size(currentframe,2)/2);
yshift = y-round(size(currentframe,1)/2);

if length(sisterList) >= currentFrame && ~isempty(sisterList{nextframe}) && ~isempty(sisterList{1})
    sisExistInd = find(sisterList{nextframe}(:,1) ~= -1 & sisterList{nextframe}(:,1) ~= 0);
else
    sisExistInd = [];
end
PosInd = find(cellpath{currentFrame}(:,1)>0 & cellpath{currentFrame}(:,2)>0)';
for cell=PosInd
    if isempty(find(cell==sisExistInd,1))
        sisterList{nextframe}(cell,:) = sisterList{currentFrame}(cell,:);
    end
    cellpath{nextframe}(cell,:) = [cellpath{currentFrame}(cell,1)+xshift cellpath{currentFrame}(cell,2)+yshift];
end

DeathInd = find(cellpath{currentFrame}(:,1)==-2)';

for cell=DeathInd
    sisterList{nextframe}(cell,:) = sisterList{currentFrame}(cell,:);
    cellpath{nextframe}(cell,:) = cellpath{currentFrame}(cell,:);
end

if ~isempty(bg)
    bg{nextframe} = [bg{currentFrame}(:,1)+xshift bg{currentFrame}(:,2)+yshift];
end

updateLists(cellpath,sisterList,res_cellpath,res_sisterList,bg,handles,nextframe)

handles.cellpath = cellpath;
handles.sisterList = sisterList;
handles.bg = bg;
guidata(hObject, handles);
handles = guidata(hObject);
restoreCellList(handles,nextframe);
imshow(imadjust((currentframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
set(handles.edit_currentFrame,'String',num2str(nextframe));
set(handles.slider_frames,'Value',nextframe);

if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,nextframe,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;

% --- Executes on button press in pushbutton_duplicatethru.
function pushbutton_duplicatethru_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_duplicatethru (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end
set(handles.togglebutton_editable,'Value',0);
cellpath = handles.cellpath;
sisterList = handles.sisterList;
res_cellpath = handles.res_cellpath;
res_sisterList = handles.res_sisterList;

bg = handles.bg;

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));

currentFrame = str2num(get(handles.edit_currentFrame,'String'));
firstframe = str2num(get(handles.edit_firstframe,'String'));
lastframe   = str2num(get(handles.edit_lastframe,'String'));

switch handles.increment
    case 1
        endframe=lastframe-1;
    case -1
        endframe=firstframe+1;
end



if currentFrame~=endframe+handles.increment
    
    handles.greenflag = 1;
    guidata(hObject, handles);
    handles = guidata(hObject);
    FFTiv=handles.FFTiv;
    FFTrv=handles.FFTrv;
    IFFTiv=handles.IFFTiv;
    previousframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],currentFrame,handles.channelnames,handles.SourceF);
    previousframe = previousframe(300:end-300,300:end-300);
    currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],currentFrame+handles.increment,handles.channelnames,handles.SourceF);
    opt = detbestlength2(FFTrv,FFTiv,IFFTiv,size(currentframe),size(previousframe),isreal(currentframe),isreal(previousframe));
    
    for t = currentFrame:handles.increment:endframe
        tp=t;
        previousframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
        previousframe = previousframe(300:end-300,300:end-300);
        tp=t+handles.increment;
        currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
        imshow(imadjust((currentframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
        
        [x y maxVal] = corrMatching2(currentframe,previousframe,opt);
        
        xshift = x-round(size(currentframe,2)/2);
        yshift = y-round(size(currentframe,1)/2);
        
        
        if  length(sisterList) >= tp && ~isempty(sisterList{tp}) && ~isempty(sisterList{1})
            sisExistInd = find(sisterList{tp}(:,1) ~= -1 & sisterList{tp}(:,1) ~= 0);
        else
            sisExistInd = [];
        end
        PosInd = find(cellpath{tp-handles.increment}(:,1)>0 & cellpath{tp-handles.increment}(:,2)>0)';
        for cell=PosInd
            if isempty(find(cell==sisExistInd,1))
                sisterList{tp}(cell,:) = sisterList{tp-handles.increment}(cell,:);
            end
            cellpath{tp}(cell,:) = [cellpath{tp-handles.increment}(cell,1)+xshift cellpath{tp-handles.increment}(cell,2)+yshift];
        end
        
        DeathInd = find(cellpath{tp-handles.increment}(:,1)==-2)';
        
        for cell=DeathInd
            sisterList{tp}(cell,:) = sisterList{tp-handles.increment}(cell,:);
            cellpath{tp}(cell,:) = cellpath{tp-handles.increment}(cell,:);
        end
        
        set(handles.edit_currentFrame,'String',num2str(tp));
        set(handles.slider_frames,'Value',tp);
        
        if ~isempty(cellpath) && ~isempty(cellpath{tp})
            set(handles.listbox_cells,'String',num2str( [(1:size(cellpath{tp},1))' cellpath{tp} sisterList{tp}] ));
        end
        
        if ~isempty(res_cellpath) && ~isempty(res_cellpath{tp})
            set(handles.listbox_Restcells,'String',num2str( [(1:size(res_cellpath{tp},1))' res_cellpath{tp} res_sisterList{tp}] ));
        else
            set(handles.listbox_Restcells,'String',[]);
        end
        
        
        if ~isempty(bg)
            %updating bg points by
            bg{tp}(:,1) = bg{tp-handles.increment}(:,1)+round(xshift);
            bg{tp}(:,2) = bg{tp-handles.increment}(:,2)+round(yshift);
            set(handles.listbox_bgs,'String',num2str( [(1:size(bg{tp},1))' bg{tp}] ));
            handles.bg = bg;
        end
        guidata(hObject, handles);
        handles = guidata(hObject);
        if get(handles.checkbox_cellmarking,'Value')
            [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
            handles.p = p;
            handles.bg_p = bg_p;
            guidata(hObject, handles);
        end
        drawnow;
        handles = guidata(hObject);
        if handles.greenflag==0
            break;
        end
    end
    handles.greenflag = 1;
    handles.cellpath = cellpath;
    handles.sisterList = sisterList;
    if tp==endframe
        set(handles.edit_commu,'String',{'Optimize points and press Record to save positions for each frame.'});
    end
    guidata(hObject, handles);
    handles = guidata(hObject);
    restoreCellList(handles,tp);
end
% --- Executes when selected object is changed in uibuttongroup_trackdirection.
function uibuttongroup_trackdirection_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup_trackdirection
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'togglebutton_backward'
        handles.increment = -1;
    case 'togglebutton_forward'
        handles.increment = 1;
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton_autolevel.
function pushbutton_autolevel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_autolevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cellpath = handles.cellpath;
bg = handles.bg;
sisterList = handles.sisterList;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
c_tp = str2num(get(handles.edit_currentFrame,'String'));
set(handles.slider_frames,'Value',str2num(get(handles.edit_currentFrame,'String')));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));
if c_tp>last_tp
    tp=last_tp;
elseif c_tp<first_tp
    tp=first_tp;
else
    tp=c_tp;
end

currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
set(handles.edit_currentFrame,'String',num2str(tp));
thres = stretchlim(currentframe);
set(handles.edit_thresMin,'String',num2str(thres(1)));
set(handles.edit_thresMax,'String',num2str(thres(2)));


imshow(imadjust((currentframe),[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;



function edit_minstamp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minstamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minstamp as text
%        str2double(get(hObject,'String')) returns contents of edit_minstamp as a double


% --- Executes during object creation, after setting all properties.
function edit_minstamp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minstamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_secstamp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_secstamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_secstamp as text
%        str2double(get(hObject,'String')) returns contents of edit_secstamp as a double


% --- Executes during object creation, after setting all properties.
function edit_secstamp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_secstamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_createmovie.
function uipanel_createmovie_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_createmovie
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton_frameno'
        handles.framestamp = 1;
    case 'radiobutton_timestamp'
        handles.framestamp = 2;
    otherwise
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton_terminate.
function pushbutton_terminate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_terminate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));

cellpath = handles.cellpath;
bg = handles.bg;
sisterList = handles.sisterList;
cellpath = handles.cellpath;

selected_cell = str2num(get(handles.edit_cellNo,'String'));


c_tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

for tp=c_tp:1:last_tp
    cellpath{tp}(selected_cell,:) = [-2 -2];
    sisterList{tp}(selected_cell,:) = sisterList{tp-1}(selected_cell,:);
end
set(handles.edit_commu,'String',['Cell#' num2str(selected_cell) ' removed from tracking']);

set(handles.togglebutton_editable,'Value',0);
currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],c_tp,handles.channelnames,handles.SourceF);
imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
[p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,c_tp,selected_cell);
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,c_tp);
handles.p = p;
handles.bg_p = bg_p;
handles.cellpath = cellpath;
handles.sisterList = sisterList;
guidata(hObject, handles);



% --- Executes on button press in pushbutton_indivShift.
function pushbutton_indivShift_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_indivShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.initialframe)
    set(handles.edit_commu,'String',['No input images.']);
    return;
end
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
cellpath = handles.cellpath;
bg = handles.bg;
sisterList = handles.sisterList;

c_tp = str2num(get(handles.edit_currentFrame,'String'));
selected_cell = str2num(get(handles.edit_cellNo,'String'));

oldCoord = cellpath{c_tp}(selected_cell,:);
set(handles.togglebutton_editable,'Value',1);
[p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,c_tp,selected_cell);
wait(p{selected_cell});
newCoord = getPosition(p{selected_cell});
shifting = round(newCoord - oldCoord);
handles.p = p;
handles.bg_p = bg_p;
guidata(hObject, handles);
handles = guidata(hObject);


if ~isempty(get(handles.edit_shiftframes,'String'))
    shiftframes = str2num(get(handles.edit_shiftframes,'String'));
    for tp=shiftframes
        cellpath{tp}(selected_cell,:) = cellpath{tp}(selected_cell,:)+shifting;
    end
    set(handles.edit_commu,'String',['X Shift:' num2str(shifting(1)) ', Y Shift:' num2str(shifting(2))]);
end

set(handles.togglebutton_editable,'Value',0);
currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],c_tp,handles.channelnames,handles.SourceF);
imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
[p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,c_tp,selected_cell);
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,c_tp);
handles.p = p;
handles.bg_p = bg_p;
handles.cellpath = cellpath;
guidata(hObject, handles);





function edit_shiftframes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_shiftframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_shiftframes as text
%        str2double(get(hObject,'String')) returns contents of edit_shiftframes as a double


% --- Executes during object creation, after setting all properties.
function edit_shiftframes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_shiftframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_Restcells.
function listbox_Restcells_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_Restcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_Restcells contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_Restcells
selected_cell = get(hObject,'Value');

res_cellpath = handles.res_cellpath;
res_sisterList = handles.res_sisterList;
bg = handles.bg;

c_tp = str2num(get(handles.edit_currentFrame,'String'));

first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));
if c_tp>last_tp
    tp=last_tp;
elseif c_tp<first_tp
    tp=first_tp;
else
    tp=c_tp;
end

if res_cellpath{tp}(selected_cell,1) ~= -1 && res_cellpath{tp}(selected_cell,2) ~= -1
    
    set(handles.edit_cellNo,'String',num2str(selected_cell));
    set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',res_cellpath{tp}(selected_cell,1),res_cellpath{tp}(selected_cell,2))));
    set(handles.edit_sister,'String',num2str(res_sisterList{tp}(selected_cell,:)));
    
    cellsize = str2num(get(handles.edit_cellsize,'String'));
    row = str2num(get(handles.edit_row,'String'));
    col = str2num(get(handles.edit_col,'String'));
    field = str2num(get(handles.edit_field,'String'));
    plane = str2num(get(handles.edit_plane,'String'));
    channel= str2num(get(handles.edit_CH,'String'));
    
    currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
    set(handles.edit_currentFrame,'String',num2str(tp));
    imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
    if get(handles.checkbox_cellmarking,'Value')
        [p bg_p] = plotTrackpoints(handles,handles.cellpath,handles.sisterList,res_cellpath,res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
        handles.p = p;
        handles.bg_p = bg_p;
    end
    
    xL=max(res_cellpath{tp}(selected_cell,1)-cellsize,1);
    xR=min(res_cellpath{tp}(selected_cell,1)+cellsize,size(currentframe,2));
    yL=max(res_cellpath{tp}(selected_cell,2)-cellsize,1);
    yR=min(res_cellpath{tp}(selected_cell,2)+cellsize,size(currentframe,1));
    template = currentframe(yL:yR,xL:xR);
    
    switch handles.cellrecogmethod
        case 1
            [x1 y1 BW] = templateToCentroid(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog);
        case 2
            [x1 y1 BW] = templateToCentroid2(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog);
        case 3
            [x1 y1 BW] = templateToCentroid3(template,round(size(template,2)/2),round(size(template,1)/2),handles.maxI,handles.invertedLog,str2num(get(handles.edit_outersize,'String')));
    end
    
    edge = bwmorph(BW,'remove');
    imshow(imadjust(template),'Parent',handles.axes2);
    set(handles.axes2,'NextPlot','add');
    [y x]=find(edge);
    plot(handles.axes2,x,y,'.r','MarkerSize',1);
    plot(handles.axes2,res_cellpath{tp}(selected_cell,1)-xL,res_cellpath{tp}(selected_cell,2)-yL,'xr');
    set(handles.axes2,'NextPlot','replace');
    guidata(hObject, handles);
    drawnow;
    
    set(handles.edit_commu,'String',['Showing Cell#' num2str(selected_cell)]);
else
    set(handles.edit_commu,'String',['Cell#' num2str(selected_cell) ' does not exist at current time point.']);
    set(handles.edit_cellNo,'String',num2str(selected_cell));
    set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',res_cellpath{tp}(selected_cell,1),res_cellpath{tp}(selected_cell,2))));
    set(handles.edit_sister,'String',num2str(res_sisterList{tp}(selected_cell,:)));
    
end

% --- Executes during object creation, after setting all properties.
function listbox_Restcells_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_Restcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_reserveall.
function pushbutton_reserveall_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reserveall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cellpath = handles.cellpath;
sisterList = handles.sisterList;
res_cellpath = handles.res_cellpath;
res_sisterList = handles.res_sisterList;
bg = handles.bg;

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

if ~isempty(cellpath)
    if ~isempty(res_cellpath)
        [sel_cellpath,sel_sisterList] = removeSister(cellpath,sisterList,first_tp,last_tp,1:size(cellpath{last_tp},1));
        
        for t = first_tp:last_tp
            if ~isempty(res_cellpath)
                new_res_cellpath{t} = [res_cellpath{t};sel_cellpath{t}];
                new_res_sisterList{t} = [res_sisterList{t};sel_sisterList{t}];
            else
                new_res_cellpath{t} = sel_cellpath{t};
                new_res_sisterList{t} = sel_sisterList{t};
            end
        end
    else
        new_res_cellpath = cellpath;
        new_res_sisterList = sisterList;
    end
    new_cellpath = [];
    new_sisterList = [];
    handles.cellpath = new_cellpath;
    handles.sisterList = new_sisterList;
    handles.res_cellpath = new_res_cellpath;
    handles.res_sisterList = new_res_sisterList;
    guidata(hObject, handles);
    
    set(handles.listbox_cells,'Value',1);
    set(handles.listbox_Restcells,'Value',1);
    
    currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
    imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
    if get(handles.checkbox_cellmarking,'Value')
        [p bg_p] = plotTrackpoints(handles,new_cellpath,new_sisterList,new_res_cellpath,new_res_sisterList,bg,tp,1);
    end
    drawnow;
    updateLists(new_cellpath,new_sisterList,new_res_cellpath,new_res_sisterList,handles.bg,handles,tp);
    guidata(hObject, handles);
end




% --- Executes on button press in pushbutton_restoreindv.
function pushbutton_restoreindv_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_restoreindv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cellpath = handles.cellpath;
sisterList = handles.sisterList;
res_cellpath = handles.res_cellpath;
res_sisterList = handles.res_sisterList;
bg = handles.bg;

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

selected_cell = str2num(get(handles.edit_cellNo,'String'));

oriSis = res_sisterList{last_tp}(selected_cell,:);
noSisInd = find(oriSis==-1,1,'first');
if ~isempty(noSisInd)
    switch noSisInd
        case 1
            moveLists = [selected_cell];
        case 2
            moveLists = [selected_cell oriSis(1)];
        case 3
            moveLists = [selected_cell oriSis(1:2)];
    end
else
    moveLists = [selected_cell oriSis];
end
if ~isempty(res_cellpath)
    % when the selected cells has no sisters
    if isempty(setdiff(moveLists,selected_cell))
        
        for t = first_tp:last_tp
            %get rid of selected cells
            restListsize = size(res_cellpath{t},1);

            belowL_rescellpath = res_cellpath{t}(1:(selected_cell-1),:);
            belowL_ressisterList = res_sisterList{t}(1:(selected_cell-1),:);

            if selected_cell<restListsize
                aboveL_rescellpath = res_cellpath{t}((selected_cell+1):restListsize,:);
                aboveL_ressisterList = res_sisterList{t}((selected_cell+1):restListsize,:);
            else
                aboveL_rescellpath = [];
                aboveL_ressisterList = [];
            end
            new_res_cellpath{t} = [belowL_rescellpath;aboveL_rescellpath];
            new_res_sisterList{t} = [belowL_ressisterList;aboveL_ressisterList];
            
            if ~isempty(cellpath)
                cellListsize = size(cellpath{t},1);
                new_cellpath{t} = cellpath{t};
                new_sisterList{t} = sisterList{t};
            else
                cellListsize = 0;
            end
            
            new_cellpath{t}(cellListsize+1,:)   =  res_cellpath{t}(selected_cell,:);
            new_sisterList{t}(cellListsize+1,:) =  res_sisterList{t}(selected_cell,:);
            
        end
        cellpath = new_cellpath;
        sisterList = new_sisterList;
        
        res_cellpath = new_res_cellpath;
        res_sisterList = new_res_sisterList;

    else
        
        firstSis = moveLists;
        secondSis = moveLists;
        for s = 1:length(firstSis)
            secondSis = [secondSis setdiff(res_sisterList{last_tp}(firstSis(s),:),-1)];
        end
        thirdSis = secondSis;
        for s = 1:length(secondSis)
            thirdSis = [thirdSis setdiff(res_sisterList{last_tp}(secondSis(s),:),-1)];
        end
        moveLists = unique(thirdSis);

        [new_res_cellpath, new_res_sisterList] = recreateAllList(res_cellpath,res_sisterList,first_tp,last_tp,moveLists,0);
        [sel_cellpath,sel_sisterList] = removeSister(res_cellpath,res_sisterList,first_tp,last_tp,moveLists);
        
        
        for t = first_tp:last_tp
            if ~isempty(cellpath)
                new_cellpath{t} = [cellpath{t};sel_cellpath{t}];
                new_sisterList{t} = [sisterList{t};sel_sisterList{t}];
            else
                new_cellpath{t} = sel_cellpath{t};
                new_sisterList{t} = sel_sisterList{t};
            end

        end
    end

    handles.cellpath = new_cellpath;
    handles.sisterList = new_sisterList;
    handles.res_cellpath = new_res_cellpath;
    handles.res_sisterList = new_res_sisterList;
    guidata(hObject, handles);
    
    set(handles.listbox_cells,'Value',1);
    set(handles.listbox_Restcells,'Value',1);
    
    currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
    imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
    if get(handles.checkbox_cellmarking,'Value')
        [p bg_p] = plotTrackpoints(handles,new_cellpath,new_sisterList,new_res_cellpath,new_res_sisterList,bg,tp,1);
    end
    drawnow;
    updateLists(new_cellpath,new_sisterList,new_res_cellpath,new_res_sisterList,handles.bg,handles,tp);
    guidata(hObject, handles);
end




% --- Executes on button press in pushbutton_restoreall.
function pushbutton_restoreall_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_restoreall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cellpath = handles.cellpath;
sisterList = handles.sisterList;
res_cellpath = handles.res_cellpath;
res_sisterList = handles.res_sisterList;
bg = handles.bg;

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

if ~isempty(res_cellpath)
    if ~isempty(cellpath) 
        [sel_cellpath,sel_sisterList] = removeSister(res_cellpath,res_sisterList,first_tp,last_tp,1:size(res_cellpath{last_tp},1));
        for t = first_tp:last_tp
            if ~isempty(cellpath)
                new_cellpath{t} = [cellpath{t};sel_cellpath{t}];
                new_sisterList{t} = [sisterList{t};sel_sisterList{t}];
            else
                new_cellpath{t} = sel_cellpath{t};
                new_sisterList{t} = sel_sisterList{t};
            end
        end
        
    else
        new_cellpath = res_cellpath;
        new_sisterList = res_sisterList;
    end
    new_res_cellpath = [];
    new_res_sisterList = [];
    handles.cellpath = new_cellpath;
    handles.sisterList = new_sisterList;
    handles.res_cellpath = new_res_cellpath;
    handles.res_sisterList = new_res_sisterList;
    guidata(hObject, handles);
    
    set(handles.listbox_cells,'Value',1);
    set(handles.listbox_Restcells,'Value',1);
    
    currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
    imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
    if get(handles.checkbox_cellmarking,'Value')
        [p bg_p] = plotTrackpoints(handles,new_cellpath,new_sisterList,new_res_cellpath,new_res_sisterList,bg,tp,1);
    end
    drawnow;
    updateLists(new_cellpath,new_sisterList,new_res_cellpath,new_res_sisterList,handles.bg,handles,tp);
    guidata(hObject, handles);
end




% --- Executes on button press in pushbutton_cleanList.
function pushbutton_cleanList_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cleanList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

set(handles.edit_cellNo,'String','1');
set(handles.listbox_cells,'Value',1);
set(handles.listbox_Restcells,'Value',1);

cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg=handles.bg;

if ~isempty(cellpath)
    [new_cellpath, new_sisterList] = recreateAllList(cellpath,sisterList,first_tp,last_tp,[],1);
    handles.cellpath   = new_cellpath;
    handles.sisterList = new_sisterList;
    
    guidata(hObject, handles);
    handles = guidata(hObject);
    updateLists(new_cellpath,new_sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,tp);
end




set(handles.edit_commu,'String','Done cleaning list.');
currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,new_cellpath,new_sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;


function [new_cellpath, new_sisterList] = recreateAllList(cellpath,sisterList,first_tp,last_tp,excludeInd,marriageLog)

sis_cellpath = cell(last_tp,1);
sis_sisterList = cell(last_tp,1);
% Determine cells with sisters
withSisInd = setdiff(find(sisterList{last_tp}(:,1)~=-1),excludeInd);

if ~isempty(withSisInd)
    
    cInd=1;
    loopInd=1;
    lastInd(loopInd)=1;
    while ~isempty(withSisInd)
        firstSis = withSisInd(1)
        secondSis = sisterList{last_tp}(withSisInd(1),1);
        sis_gInd = find(sisterList{last_tp}(:,1)==firstSis | sisterList{last_tp}(:,1)==secondSis);
        
        for t = last_tp:-1:first_tp
            if  t~=last_tp
                sis_sisterList{t} = -1*ones(size(sis_sisterList{last_tp}));
            else
                
                cInd=lastInd(loopInd);
                
                for s=1:length(sis_gInd)
                    sis_table{sis_gInd(s)} = s+cInd-1;
                end
                
                for s=1:length(sis_gInd)
                    sis_cellpath{t}(cInd,:)   = cellpath{t}(sis_gInd(s),:);
                    oldList = sisterList{t}(sis_gInd(s),:);
                    posL = find(oldList ~= -1);
                    
                    if ~isempty(posL)
                        switch length(posL)
                            case 1
                                sis_sisterList{t}(cInd,:) = [sis_table{oldList(1)} -1 -1];
                            case 2
                                sis_sisterList{t}(cInd,:) = [sis_table{oldList(1)} sis_table{oldList(2)} -1];
                            case 3
                                sis_sisterList{t}(cInd,:) = [sis_table{oldList(1)} sis_table{oldList(2)} sis_table{oldList(3)}];
                        end
                    end
                    cInd = cInd+1;
                end
                
            end
            
        end
        loopInd=loopInd+1;
        lastInd(loopInd) = cInd;
        withSisInd = setdiff(withSisInd,sis_gInd);
        clear sis_table;
        
    end
    
    cInd=1;
    loopInd=1;
    lastInd(loopInd)=1;
    withSisInd = setdiff(find(sisterList{last_tp}(:,1)~=-1),excludeInd);
    while ~isempty(withSisInd)
        firstSis = withSisInd(1);
        secondSis = sisterList{last_tp}(withSisInd(1),1);
        sis_gInd = find(sisterList{last_tp}(:,1)==firstSis | sisterList{last_tp}(:,1)==secondSis);
        
        for t = last_tp:-1:first_tp
            
            cInd=lastInd(loopInd);
            
            for s=1:length(sis_gInd)
                sis_table{sis_gInd(s)} = s+cInd-1;
            end
            
            for s=1:length(sis_gInd)
                sis_cellpath{t}(cInd,:)   = cellpath{t}(sis_gInd(s),:);
                oldList = sisterList{t}(sis_gInd(s),:);
                posL = find(oldList ~= -1);
                if ~isempty(posL)
                    switch length(posL)
                        case 1
                            sis_sisterList{t}(cInd,:) = [sis_table{oldList(1)} -1 -1];
                        case 2
                            sis_sisterList{t}(cInd,:) = [sis_table{oldList(1)} sis_table{oldList(2)} -1];
                        case 3
                            sis_sisterList{t}(cInd,:) = [sis_table{oldList(1)} sis_table{oldList(2)} sis_table{oldList(3)}];
                    end
                end
                cInd = cInd+1;
            end
            
        end
        loopInd=loopInd+1;
        lastInd(loopInd) = cInd;
        withSisInd = setdiff(withSisInd,sis_gInd);
        clear sis_table;
    end
    
end

% Determine cells without sisters
noSisInd = setdiff(find(sisterList{last_tp}(:,1)==-1 & cellpath{last_tp}(:,1)~=-1),excludeInd);

for t = first_tp:last_tp
    nosis_cellpath{t}   = cellpath{t}(noSisInd,:);
    nosis_sisterList{t} = sisterList{t}(noSisInd,:) ;
end

if marriageLog
    [single_cellpath single_sisterList couple_cellpath couple_sisterList] = massWedding(nosis_cellpath,nosis_sisterList,first_tp,last_tp);
    if ~isempty(couple_cellpath)
        for t = first_tp:last_tp
            if ~isempty(sis_cellpath)
                oldsissize = size(sis_cellpath{t},1);
                new_cellpath{t} =   sis_cellpath{t};
                new_sisterList{t} = sis_sisterList{t};
            else
                oldsissize = 0;
            end
            
            addsissize = size(couple_cellpath{t},1);
            for c = 1 : addsissize
                new_cellpath{t}(oldsissize+c,:)   =  couple_cellpath{t}(c,:);
                old_sisters = couple_sisterList{t}(c,:);
                PosSis = find(old_sisters>0);
                new_sisterList{t}(oldsissize+c,:) =  old_sisters;
                if ~isempty(PosSis)
                    for ss=1:length(PosSis)
                        new_sisterList{t}(oldsissize+c,PosSis(ss)) = old_sisters(PosSis(ss))+oldsissize;
                    end
                end
                
            end
        end
        sis_cellpath = new_cellpath;
        sis_sisterList = new_sisterList;
    end

    for t = first_tp:last_tp
        new_cellpath{t} =[sis_cellpath{t};single_cellpath{t}];
        new_sisterList{t} = [sis_sisterList{t};single_sisterList{t}];
    end
    
else
    for t = first_tp:last_tp
        new_cellpath{t} =   [sis_cellpath{t};nosis_cellpath{t}];
        new_sisterList{t} = [sis_sisterList{t};nosis_sisterList{t}];
    end
end


function [single_cellpath single_sisterList couple_cellpath couple_sisterList] = massWedding(nosis_cellpath,nosis_sisterList,first_tp,last_tp)
single_cellpath   = nosis_cellpath;
single_sisterList = nosis_sisterList;

couple_cellpath   = cell(last_tp,1);
couple_sisterList = cell(last_tp,1);

Rad = 2;
%create distance matrix for first time point
D = pdist(nosis_cellpath{first_tp},'euclidean');
disMat = squareform(D);

noSisList = [];
SisList = [];
out_Ind = 1;
leftInd = 1:size(disMat,1)';
while ~isempty(leftInd)
    sis1Ind = find(disMat(leftInd(1),:)<Rad);
    sisInd = find(sis1Ind~=leftInd(1));
    if isempty(sisInd)
        noSisList = [noSisList;leftInd(1)];
        cList = [leftInd(1)];
    else
        if length(sisInd) == 1
            ind_dist = zeros(1,last_tp);
            for t=first_tp:last_tp
                ind_dist(t) = pdist([nosis_cellpath{t}(leftInd(1),:);nosis_cellpath{t}(sis1Ind(sisInd(1)),:)]);
            end
            diff_dist = diff(ind_dist);
            jumpInd = find(diff_dist>0,1,'first');
            if isempty(jumpInd)
                splitF = 1;
            else
                splitF = jumpInd+1;
            end
            for t=first_tp:splitF-1
                couple_sisterList{t}(out_Ind,:) = [-1 -1 -1];
                couple_sisterList{t}(out_Ind+1,:) = [-1 -1 -1];
                couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(leftInd(1),:);
                couple_cellpath{t}(out_Ind+1,:) = [-1 -1];
            end
            for t=splitF:last_tp
                couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                couple_sisterList{t}(out_Ind+1,:) = [out_Ind -1 -1];
                couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(leftInd(1),:);
                couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(sis1Ind(sisInd(1)),:);
            end
            SisList = [SisList;leftInd(1);sis1Ind(sisInd(1))];
            cList = [leftInd(1);sis1Ind(sisInd(1))];
            out_Ind=out_Ind+2;
            
        elseif length(sisInd) == 2
            cellNos = [leftInd(1) sis1Ind(sisInd)]
            ind_dist = zeros(1,last_tp);
            c_disMat = zeros(length(sisInd)+1,length(sisInd)+1,last_tp-1);
            for t=first_tp:(last_tp-1)
                c_mat1 = [nosis_cellpath{t}(leftInd(1),:)];
                for s = 1:length(sisInd)
                    c_mat1 = [c_mat1;nosis_cellpath{t}(sis1Ind(sisInd(s)),:)];
                end
                c_mat2 = [nosis_cellpath{t+1}(leftInd(1),:)];
                for s = 1:length(sisInd)
                    c_mat2 = [c_mat2;nosis_cellpath{t+1}(sis1Ind(sisInd(s)),:)];
                end
                c_disMat(:,:,t) = squareform(pdist(c_mat2))-squareform(pdist(c_mat1));
                clear c_mat;
            end
            TimeDistance = zeros(size(c_disMat,1),size(c_disMat,2));
            for i=1:size(c_disMat,1)
                for j=1:size(c_disMat,2)
                    if i~=j
                        TimeDistance(i,j) = find(squeeze(c_disMat(i,j,:)),1,'first');
                    end
                end
            end
            DivTime = unique(TimeDistance(TimeDistance~=0));
            
            switch length(DivTime)
                case 1               
                    for t=first_tp:DivTime
                        couple_sisterList{t}(out_Ind,:) = [-1 -1 -1];
                        couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(cellNos(1),:);
                        for i=2:length(cellNos)
                            couple_sisterList{t}(out_Ind+i-1,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+i-1,:) = [-1 -1];
                        end
                        
                    end
                    for t=DivTime+1:last_tp
                        couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                        couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(cellNos(1),:);
                        
                        for i=2:length(cellNos)
                            couple_sisterList{t}(out_Ind+i-1,:) = [out_Ind -1 -1];
                            couple_cellpath{t}(out_Ind+i-1,:) = nosis_cellpath{t}(cellNos(i),:);
                        end
                    end
                    
                case 2
                    DivTime1 = min(DivTime);
                    DivTime2 = max(DivTime);
                    [r,c] = find(TimeDistance==DivTime1);
                    TCounts = histc(r,1:size(c_disMat,1));
                    Div1Ind = find(TCounts==max(TCounts));
                    Div1Cell = cellNos(Div1Ind);
                    Div2Cells = setdiff(cellNos,Div1Cell);
                    
                    for t=first_tp:(DivTime1)
                        couple_sisterList{t}(out_Ind,:) = [-1 -1 -1];
                        couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(Div1Cell,:);
                        couple_sisterList{t}(out_Ind+1,:) = [-1 -1 -1];
                        couple_cellpath{t}(out_Ind+1,:) = [-1 -1];
                        couple_sisterList{t}(out_Ind+2,:) = [-1 -1 -1];
                        couple_cellpath{t}(out_Ind+2,:) = [-1 -1];
                        
                    end
                    for t=(DivTime1+1):(DivTime2)
                        couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                        couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(Div1Cell,:);
                        couple_sisterList{t}(out_Ind+1,:) = [out_Ind -1 -1];
                        couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(Div2Cells(1),:);
                        couple_sisterList{t}(out_Ind+2,:) = [out_Ind -1 -1];
                        couple_cellpath{t}(out_Ind+2,:) = [-1 -1];
                    end
                    
                    for t = (DivTime2+1):last_tp
                        couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                        couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(Div1Cell,:);
                        couple_sisterList{t}(out_Ind+1,:) = [out_Ind out_Ind+2 -1];
                        couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(Div2Cells(1),:);
                        couple_sisterList{t}(out_Ind+2,:) = [out_Ind out_Ind+1 -1];
                        couple_cellpath{t}(out_Ind+2,:) = nosis_cellpath{t}(Div2Cells(2),:);
                    end
                    
            end
            out_Ind=out_Ind+3;
            cList = [leftInd(1);sis1Ind(sisInd(1));sis1Ind(sisInd(2))];
        elseif length(sisInd) == 3
            cellNos = [leftInd(1) sis1Ind(sisInd)]
            ind_dist = zeros(1,last_tp);
            c_disMat = zeros(length(sisInd)+1,length(sisInd)+1,last_tp-1);
            for t=first_tp:(last_tp-1)
                c_mat1 = [nosis_cellpath{t}(leftInd(1),:)];
                for s = 1:length(sisInd)
                    c_mat1 = [c_mat1;nosis_cellpath{t}(sis1Ind(sisInd(s)),:)];
                end
                c_mat2 = [nosis_cellpath{t+1}(leftInd(1),:)];
                for s = 1:length(sisInd)
                    c_mat2 = [c_mat2;nosis_cellpath{t+1}(sis1Ind(sisInd(s)),:)];
                end
                c_disMat(:,:,t) = squareform(pdist(c_mat2))-squareform(pdist(c_mat1));
                clear c_mat;
            end
            TimeDistance = zeros(size(c_disMat,1),size(c_disMat,2));
            for i=1:size(c_disMat,1)
                for j=1:size(c_disMat,2)
                    if i~=j
                        TimeDistance(i,j) = find(squeeze(c_disMat(i,j,:)),1,'first');
                    end
                end
            end
            DivTime = unique(TimeDistance(TimeDistance~=0));
            
            switch length(DivTime)
                case 1
                    for t=first_tp:DivTime
                        couple_sisterList{t}(out_Ind,:) = [-1 -1 -1];
                        couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(cellNos(1),:);
                        for i=2:length(cellNos)
                            couple_sisterList{t}(out_Ind+i-1,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+i-1,:) = [-1 -1];
                        end
                    end
                    for t=DivTime+1:last_tp
                        couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                        couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(cellNos(1),:);
                        
                        for i=2:length(cellNos)
                            couple_sisterList{t}(out_Ind+i-1,:) = [out_Ind -1 -1];
                            couple_cellpath{t}(out_Ind+i-1,:) = nosis_cellpath{t}(cellNos(i),:);
                        end
                    end
                    
                case 2
                    DivTime1 = min(DivTime);
                    DivTime2 = max(DivTime);
                    [r,c] = find(TimeDistance==DivTime1);
                    TCounts = histc(r,1:size(c_disMat,1));
                    Div1Ind = find(TCounts==max(TCounts));
                    Div1Cell = cellNos(Div1Ind);
                    Div2Cells = setdiff(cellNos,Div1Cell);
                    
                    for t=first_tp:(DivTime1)
                        couple_sisterList{t}(out_Ind,:) = [-1 -1 -1];
                        couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(Div1Cell,:);
                        couple_sisterList{t}(out_Ind+1,:) = [-1 -1 -1];
                        couple_cellpath{t}(out_Ind+1,:) = [-1 -1];
                        couple_sisterList{t}(out_Ind+2,:) = [-1 -1 -1];
                        couple_cellpath{t}(out_Ind+2,:) = [-1 -1];
                    end
                    for t=(DivTime1+1):(DivTime2)
                        couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                        couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(Div1Cell,:);
                        couple_sisterList{t}(out_Ind+1,:) = [out_Ind -1 -1];
                        couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(Div2Cells(1),:);
                        couple_sisterList{t}(out_Ind+2,:) = [out_Ind -1 -1];
                        couple_cellpath{t}(out_Ind+2,:) = [-1 -1];
                    end
                    for t=(DivTime2+1):last_tp
                        couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                        couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(Div1Cell,:);
                        couple_sisterList{t}(out_Ind+1,:) = [out_Ind out_Ind+2 -1];
                        couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(Div2Cells(1),:);
                        couple_sisterList{t}(out_Ind+2,:) = [out_Ind out_Ind+1 -1];
                        couple_cellpath{t}(out_Ind+2,:) = nosis_cellpath{t}(Div2Cells(2),:);
                    end
                    
                case 3
                    DivTime1 = min(DivTime);
                    DivTime3 = max(DivTime);
                    DivTime2 = setdiff(DivTime,[DivTime1 DivTime3]);
                    [r,c] = find(TimeDistance==DivTime1);
                    TCounts = histc(r,1:size(c_disMat,1));
                    
                    Div1Ind = find(TCounts==max(TCounts));
                    
                    if length(Div1Ind) ~= 1
                        
                        [r,c] = find(TimeDistance==DivTime2);
                        TCounts = histc(r,1:size(c_disMat,1));
                        Div2Ind = find(TCounts==max(TCounts));
                        Div2Cells = cellNos(Div2Ind);
                        
                        [r,c] = find(TimeDistance==DivTime3);
                        TCounts = histc(r,1:size(c_disMat,1));
                        Div3Ind = find(TCounts==max(TCounts));
                        Div3Cells = cellNos(Div3Ind);
                        
                        ParentCell  = Div2Cells(1);
                        Sister1Cell = Div3Cells(1);
                        Sister2Cell = Div2Cells(2);
                        Sister3Cell = Div3Cells(2);
                        
                        for t=first_tp:DivTime1
                            couple_sisterList{t}(out_Ind,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(ParentCell,:);
                            couple_sisterList{t}(out_Ind+1,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+1,:) = [-1 -1];
                            couple_sisterList{t}(out_Ind+2,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+2,:) = [-1 -1];
                            couple_sisterList{t}(out_Ind+3,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+3,:) = [-1 -1];
                        end
                        for t=(DivTime1+1):DivTime2
                            couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                            couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(ParentCell,:);
                            couple_sisterList{t}(out_Ind+1,:) = [out_Ind -1 -1];
                            couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(Sister1Cell,:);
                            couple_sisterList{t}(out_Ind+2,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+2,:) = [-1 -1];
                            couple_sisterList{t}(out_Ind+3,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+3,:) = [-1 -1];
                        end
                        for t=(DivTime2+1):DivTime3
                            couple_sisterList{t}(out_Ind,:) = [out_Ind+1 out_Ind+2 -1];
                            couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(ParentCell,:);
                            couple_sisterList{t}(out_Ind+1,:) = [out_Ind -1 -1];
                            couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(Sister1Cell,:);
                            couple_sisterList{t}(out_Ind+2,:) = [out_Ind+1 out_Ind -1];
                            couple_cellpath{t}(out_Ind+2,:) = nosis_cellpath{t}(Sister2Cell,:);
                            couple_sisterList{t}(out_Ind+3,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+3,:) = [-1 -1];
                        end
                        for t=(DivTime3+1):last_tp
                            couple_sisterList{t}(out_Ind,:) = [out_Ind+1 out_Ind+2 -1];
                            couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(ParentCell,:);
                            couple_sisterList{t}(out_Ind+1,:) = [out_Ind out_Ind+3 -1];
                            couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(Sister1Cell,:);
                            couple_sisterList{t}(out_Ind+2,:) = [out_Ind+1 out_Ind -1];
                            couple_cellpath{t}(out_Ind+2,:) = nosis_cellpath{t}(Sister2Cell,:);
                            couple_sisterList{t}(out_Ind+3,:) = [out_Ind out_Ind+1 -1];
                            couple_cellpath{t}(out_Ind+3,:) = nosis_cellpath{t}(Sister3Cell,:);
                        end
                        
                    else
                        ParentCell = cellNos(Div1Ind);
                        
                        [r,c] = find(TimeDistance==DivTime2);
                        TCounts = histc(r,1:size(c_disMat,1));
                        Div2Ind = find(TCounts==max(TCounts));
                        Sister1Cell = cellNos(Div2Ind);
                        
                        [r,c] = find(TimeDistance==DivTime3);
                        TCounts = histc(r,1:size(c_disMat,1));
                        Div3Ind = find(TCounts==max(TCounts));
                        Div3Cells = cellNos(Div3Ind);
                        
                        Sister2Cell = Div3Cells(1);
                        Sister3Cell = Div3Cells(2);
                        
                        
                        for t=first_tp:DivTime1
                            couple_sisterList{t}(out_Ind,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(ParentCell,:);
                            couple_sisterList{t}(out_Ind+1,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+1,:) = [-1 -1];
                            couple_sisterList{t}(out_Ind+2,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+2,:) = [-1 -1];
                            couple_sisterList{t}(out_Ind+3,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+3,:) = [-1 -1];
                        end
                        for t=(DivTime1+1):DivTime2
                            couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                            couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(ParentCell,:);
                            couple_sisterList{t}(out_Ind+1,:) = [out_Ind -1 -1];
                            couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(Sister1Cell,:);
                            couple_sisterList{t}(out_Ind+2,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+2,:) = [-1 -1];
                            couple_sisterList{t}(out_Ind+3,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+3,:) = [-1 -1];
                        end
                        for t=(DivTime2+1):DivTime3
                            couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                            couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(ParentCell,:);
                            couple_sisterList{t}(out_Ind+1,:) = [out_Ind out_Ind+2 -1];
                            couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(Sister1Cell,:);
                            couple_sisterList{t}(out_Ind+2,:) = [out_Ind out_Ind+1 -1];
                            couple_cellpath{t}(out_Ind+2,:) = nosis_cellpath{t}(Sister2Cell,:);
                            couple_sisterList{t}(out_Ind+3,:) = [-1 -1 -1];
                            couple_cellpath{t}(out_Ind+3,:) = [-1 -1];
                        end
                        for t=(DivTime3+1):last_tp
                            couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                            couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(ParentCell,:);
                            couple_sisterList{t}(out_Ind+1,:) = [out_Ind out_Ind+2 -1];
                            couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(Sister1Cell,:);
                            couple_sisterList{t}(out_Ind+2,:) = [out_Ind out_Ind+1 out_Ind+3];
                            couple_cellpath{t}(out_Ind+2,:) = nosis_cellpath{t}(Sister2Cell,:);
                            couple_sisterList{t}(out_Ind+3,:) = [out_Ind out_Ind+1 out_Ind+2];
                            couple_cellpath{t}(out_Ind+3,:) = nosis_cellpath{t}(Sister3Cell,:);
                        end
                    end

            end
            out_Ind=out_Ind+4;
            cList = [leftInd(1);sis1Ind(sisInd(1));sis1Ind(sisInd(2));sis1Ind(sisInd(3))];
        end
    end
    leftInd = setdiff(leftInd,cList);
    
end

for t=first_tp:last_tp
    single_cellpath{t} = nosis_cellpath{t}(noSisList,:);
    single_sisterList{t} = nosis_sisterList{t}(noSisList,:);
end

% --- Executes on button press in pushbutton_makesisters.
function pushbutton_makesisters_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_makesisters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg = handles.bg;
c_tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

joinSis = str2num(get(handles.edit_joiningsisters,'String'));

main_cell = joinSis(1);
oriSis = sisterList{last_tp}(main_cell,:);
noSisInd = find(oriSis==-1,1,'first');
if ~isempty(noSisInd)
    switch noSisInd
        case 1
            main_cell_sis = [];
        case 2
            main_cell_sis = oriSis(1);
        case 3
            main_cell_sis = oriSis(1:2);
    end
else
    set(handles.edit_commu,'String',['Sister list is already full for ' num2str(main_cell)]);
    return;
end

sec_cell = joinSis(2);
oriSis = sisterList{last_tp}(sec_cell,:);
noSisInd = find(oriSis==-1,1,'first');
if ~isempty(noSisInd)
    switch noSisInd
        case 1
            sec_cell_sis = [];
        case 2
            sec_cell_sis = oriSis(1);
        case 3
            sec_cell_sis = oriSis(1:2);
    end
else
    set(handles.edit_commu,'String',['Sister list is already full for ' num2str(sec_cell)]);
    return;
end
if isempty(main_cell_sis) && isempty(sec_cell_sis)
    for t = first_tp:c_tp-1
        sisterList{t}(main_cell,:) = [-1 -1 -1];
        
        for cell = joinSis(2:end)
            cellpath{t}(cell,:)   = [-1 -1];
            sisterList{t}(cell,:) = [-1 -1 -1];
        end
    end
    for t = c_tp:last_tp
        sisterList{t}(main_cell,:) = [joinSis(2) -1 -1];
        
        for cell = joinSis(2:end)
            sisterList{t}(cell,:) = [main_cell -1 -1];
        end
    end
elseif isempty(main_cell_sis) && ~isempty(sec_cell_sis)
    
    firstSis = sec_cell_sis;
    secondSis = sec_cell_sis;
    for s = 1:length(firstSis)
        secondSis = [secondSis setdiff(sisterList{last_tp}(firstSis(s),:),-1)];
    end
    thirdSis = secondSis;
    for s = 1:length(secondSis)
        thirdSis = [thirdSis setdiff(sisterList{last_tp}(secondSis(s),:),-1)];
    end
    sec_SisList = unique(thirdSis);
    
    
    for t = first_tp:c_tp-1
        sisterList{t}(main_cell,:) = [-1 -1 -1];
        cellpath{t}(sec_cell,:)    = [-1 -1];
        sisterList{t}(sec_cell,:)  = [-1 -1 -1];
    end
    for t = c_tp:last_tp
        sisterList{t}(main_cell,:) = [sec_cell -1 -1];
        for s = 1:length(sec_SisList)
            oriSis = sisterList{t}(sec_SisList(s),:);
            noSisInd = find(oriSis==-1,1,'first');
            if ~isempty(noSisInd)
                switch noSisInd
                    case 1
                        sisterList{t}(sec_SisList(s),:) = [main_cell -1 -1];
                    case 2
                        sisterList{t}(sec_SisList(s),:) = [main_cell oriSis(1) -1];
                    case 3
                        sisterList{t}(sec_SisList(s),:) = [main_cell oriSis(1:2) ];
                    otherwise
                        set(handles.edit_commu,'String',['Sister list is already full for ' num2str(sec_SisList(s))]);
                        return;
                end
            end
        end
    end
elseif ~isempty(main_cell_sis) && isempty(sec_cell_sis)
    firstSis = main_cell_sis;
    secondSis = main_cell_sis;
    for s = 1:length(firstSis)
        secondSis = [secondSis setdiff(sisterList{last_tp}(firstSis(s),:),-1)];
    end
    thirdSis = secondSis;
    for s = 1:length(secondSis)
        thirdSis = [thirdSis setdiff(sisterList{last_tp}(secondSis(s),:),-1)];
    end
    main_SisList = unique(thirdSis);

    for t = first_tp:c_tp-1
        sisterList{t}(main_cell,:) = [-1 -1 -1];
        cellpath{t}(sec_cell,:)    = [-1 -1];
        sisterList{t}(sec_cell,:)  = [-1 -1 -1];
    end
    
    for t = c_tp:last_tp
        sisterList{t}(sec_cell,:) = [main_cell -1 -1];
        for s = 1:length(main_SisList)
            oriSis = sisterList{t}(main_SisList(s),:);
            noSisInd = find(oriSis==-1,1,'first');
            if ~isempty(noSisInd)
                switch noSisInd
                    case 1
                        sisterList{t}(main_SisList(s),:) = [sec_cell -1 -1];
                    case 2
                        sisterList{t}(main_SisList(s),:) = [sec_cell oriSis(1) -1];
                    case 3
                        sisterList{t}(main_SisList(s),:) = [sec_cell oriSis(1:2) ];
                    otherwise
                        set(handles.edit_commu,'String',['Sister list is already full for ' num2str(main_SisList(s))]);
                        return;
                end
            end
        end
    end
        
else
    firstSis = main_cell_sis;
    secondSis = main_cell_sis;
    for s = 1:length(firstSis)
        secondSis = [secondSis setdiff(sisterList{last_tp}(firstSis(s),:),-1)];
    end
    thirdSis = secondSis;
    for s = 1:length(secondSis)
        thirdSis = [thirdSis setdiff(sisterList{last_tp}(secondSis(s),:),-1)];
    end
    main_SisList = unique(thirdSis);
    
    firstSis = sec_cell_sis;
    secondSis = sec_cell_sis;
    for s = 1:length(firstSis)
        secondSis = [secondSis setdiff(sisterList{last_tp}(firstSis(s),:),-1)];
    end
    thirdSis = secondSis;
    for s = 1:length(secondSis)
        thirdSis = [thirdSis setdiff(sisterList{last_tp}(secondSis(s),:),-1)];
    end
    sec_SisList = unique(thirdSis);
    
    
    for t = first_tp:c_tp-1
        sisterList{t}(main_cell,:) = [-1 -1 -1];
        cellpath{t}(sec_cell,:)    = [-1 -1];
        sisterList{t}(sec_cell,:)  = [-1 -1 -1];
    end
    for t = c_tp:last_tp
        
        for s = 1:length(main_SisList)
            oriSis = sisterList{t}(main_SisList(s),:);
            noSisInd = find(oriSis==-1,1,'first');
            if ~isempty(noSisInd)
                switch noSisInd
                    case 1
                        sisterList{t}(main_SisList(s),:) = [sec_cell -1 -1];
                    case 2
                        sisterList{t}(main_SisList(s),:) = [sec_cell oriSis(1) -1];
                    case 3
                        sisterList{t}(main_SisList(s),:) = [sec_cell oriSis(1:2) ];
                    otherwise
                        set(handles.edit_commu,'String',['Sister list is already full for ' num2str(main_SisList(s))]);
                        return;
                end
            end
        end
        
        for s = 1:length(sec_SisList)
            oriSis = sisterList{t}(sec_SisList(s),:);
            noSisInd = find(oriSis==-1,1,'first');
            if ~isempty(noSisInd)
                switch noSisInd
                    case 1
                        sisterList{t}(sec_SisList(s),:) = [main_cell -1 -1];
                    case 2
                        sisterList{t}(sec_SisList(s),:) = [main_cell oriSis(1) -1];
                    case 3
                        sisterList{t}(sec_SisList(s),:) = [main_cell oriSis(1:2) ];
                    otherwise
                        set(handles.edit_commu,'String',['Sister list is already full for ' num2str(sec_SisList(s))]);
                        return;
                end
            end
        end
    end
end

handles.cellpath   = cellpath;
handles.sisterList = sisterList;

currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],c_tp,handles.channelnames,handles.SourceF);
imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,c_tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,c_tp);
guidata(hObject, handles);


function edit_joiningsisters_Callback(hObject, eventdata, handles)
% hObject    handle to edit_joiningsisters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_joiningsisters as text
%        str2double(get(hObject,'String')) returns contents of edit_joiningsisters as a double


% --- Executes during object creation, after setting all properties.
function edit_joiningsisters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_joiningsisters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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


% --- Executes on button press in pushbutton_clearsisters.
function pushbutton_clearsisters_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clearsisters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

cellpath = handles.cellpath;
sisterList = handles.sisterList;
bg=handles.bg;

if ~isempty(cellpath)

    [new_cellpath,new_sisterList] = removeSister(cellpath,sisterList,first_tp,last_tp,1:length(cellpath{last_tp}));
    handles.cellpath   = new_cellpath;
    handles.sisterList = new_sisterList;
    
    guidata(hObject, handles);
    handles = guidata(hObject);
    updateLists(new_cellpath,new_sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,tp);
end

set(handles.edit_commu,'String','Removed all sister information');
currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],tp,handles.channelnames,handles.SourceF);
imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
if get(handles.checkbox_cellmarking,'Value')
    [p bg_p] = plotTrackpoints(handles,handles.cellpath,handles.sisterList,handles.res_cellpath,handles.res_sisterList,bg,tp,str2num(get(handles.edit_cellNo,'String')));
    handles.p = p;
    handles.bg_p = bg_p;
    guidata(hObject, handles);
end
drawnow;



function [new_cellpath,new_sisterList] = removeSister(cellpath,sisterList,first_tp,last_tp,testList)
sis_cellpath = cell(last_tp,1);
sis_sisterList = cell(last_tp,1);
cInd=1;
% Determine cells with sisters
withSisInd = intersect(find(sisterList{last_tp}(:,1)~=-1),testList);
% Determine cells without sisters

while ~isempty(withSisInd) 
    SisList = withSisInd(1);
    firstSis = setdiff(sisterList{last_tp}(withSisInd(1),:),-1);
    secondSis = firstSis;
    for s = 1:length(firstSis)
        secondSis = [secondSis setdiff(sisterList{last_tp}(firstSis(s),:),-1)];
    end
    thirdSis = secondSis;
    for s = 1:length(secondSis)
        thirdSis = [thirdSis setdiff(sisterList{last_tp}(secondSis(s),:),-1)];
    end
    
    SisList = unique(thirdSis)
    
    for s=1:length(SisList)
        
        for t = first_tp:last_tp
            if cellpath{t}(SisList(s),1) ~= -1 && cellpath{t}(SisList(s),2) ~= -1
                sis_cellpath{t}(cInd,:)   = cellpath{t}(SisList(s),:);
                sis_sisterList{t}(cInd,:) = [-1 -1 -1];
            else
                havecoordInd = find(cellpath{t}(SisList,1) ~= -1 & cellpath{t}(SisList,2) ~= -1);
                if length(havecoordInd) == 1
                    sis_cellpath{t}(cInd,:)   = cellpath{t}(SisList(havecoordInd),:);
                    sis_sisterList{t}(cInd,:) = [-1 -1 -1];
                else
                    oriSis = sisterList{last_tp}(SisList(s),:);
                    NegOneInd = find(oriSis==-1,1,'first');
                    if ~isempty(NegOneInd)
                        mySis = oriSis(1:(NegOneInd-1));
                    else
                        mySis = oriSis;
                    end
                    for ms = length(mySis):-1:1
                        if cellpath{t}(mySis(ms),1) ~= -1 || cellpath{t}(mySis(ms),2) ~= -1
                            mytrueParent = mySis(ms);
                            break
                        end
                    end
                    sis_cellpath{t}(cInd,:)   = cellpath{t}(mytrueParent,:);
                    sis_sisterList{t}(cInd,:) = [-1 -1 -1];
                end
            end
        end
        cInd = cInd+1;
    end
    
    withSisInd = setdiff(withSisInd,SisList);
end

noSisterInd = intersect(find(sisterList{last_tp}(:,1)==-1 & cellpath{last_tp}(:,1)~=-1),testList);
for t = first_tp:last_tp
    nosis_cellpath{t}   = cellpath{t}(noSisterInd,:);
    nosis_sisterList{t} = sisterList{t}(noSisterInd,:) ;
end

% Combine data
for t = first_tp:last_tp
    new_cellpath{t} =[sis_cellpath{t};nosis_cellpath{t}];
    new_sisterList{t} = [sis_sisterList{t};nosis_sisterList{t}];
end

% --- Executes on button press in pushbutton_deleteDup.
function pushbutton_deleteDup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_deleteDup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.initialframe)
    set(handles.edit_commu,'String','No input images.');
    return;
end

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channel= str2num(get(handles.edit_CH,'String'));
c_tp = str2num(get(handles.edit_currentFrame,'String'));

cellpath = handles.cellpath;
bg = handles.bg;
sisterList = handles.sisterList;
del_cells = [];
for c = 1:size(cellpath{c_tp},1)
    if isempty(find(c==del_cells,1))
        dup_ind = intersect(find(cellpath{c_tp}(c,1)==cellpath{c_tp}(:,1)),find(cellpath{c_tp}(c,2)==cellpath{c_tp}(:,2)));
        if ~isempty(dup_ind)
            non_currentInd = find(dup_ind~=c);
            del_cells = [del_cells dup_ind(non_currentInd)'];
        end
    end
end

first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

if ~isempty(cellpath)
    if length(cellpath)<last_tp
        last_tp = length(cellpath);
    end
end

for d = del_cells
    for tp=first_tp:last_tp
        cellpath{tp}(d,:) = [-1 -1];
        sisterList{tp}(d,:) = [-1 -1 -1];
    end
end
set(handles.edit_commu,'String','Done removing duplicates.');
currentframe = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane channel],c_tp,handles.channelnames,handles.SourceF);
imshow(imadjust(currentframe,[str2num(get(handles.edit_thresMin,'String')) str2num(get(handles.edit_thresMax,'String'))],[0 1]),'Parent',handles.axes1);
[p bg_p] = plotTrackpoints(handles,cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,c_tp,1);
updateLists(cellpath,sisterList,handles.res_cellpath,handles.res_sisterList,bg,handles,c_tp);
handles.p = p;
handles.bg_p = bg_p;
handles.cellpath = cellpath;
handles.sisterList = sisterList;
guidata(hObject, handles);
drawnow;
