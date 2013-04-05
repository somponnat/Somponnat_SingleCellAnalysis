function varargout = MaskGenerating(varargin)
% MASKGENERATING MATLAB code for MaskGenerating.fig
%      MASKGENERATING, by itself, creates a new MASKGENERATING or raises the existing
%      singleton*.
%
%      H = MASKGENERATING returns the handle to a new MASKGENERATING or the handle to
%      the existing singleton*.
%
%      MASKGENERATING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MASKGENERATING.M with the given input arguments.
%
%      MASKGENERATING('Property','Value',...) creates a new MASKGENERATING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MaskGenerating_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MaskGenerating_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MaskGenerating

% Last Modified by GUIDE v2.5 09-Jan-2013 10:41:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MaskGenerating_OpeningFcn, ...
                   'gui_OutputFcn',  @MaskGenerating_OutputFcn, ...
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


% --- Executes just before MaskGenerating is made visible.
function MaskGenerating_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MaskGenerating (see VARARGIN)

% Choose default command line output for MaskGenerating
handles.output = hObject;

filetype = []; 
trackinginfo = [];  
imagelocation = []; 
channelnames = []; 
fileformat = [];   
currentcell = 1;
stageLoc = 1;

handles.currentnucmask = [];
handles.currentcytomask = [];
handles.currentcellmask = [];
handles.selected_cells = [];
handles.ind_cellpath = [];

switch length(varargin)
    case 1
        filetype = varargin{1}; 
    case 2
        filetype = varargin{1}; 
        if ischar(varargin{2}) & filetype == 3
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
        set(handles.edit_template,'String',num2str(trackinginfo(1)));
        set(handles.edit_nomin,'String',num2str(trackinginfo(2)));
        set(handles.edit_denomin,'String',num2str(trackinginfo(3)));
        set(handles.edit_firsttp,'String',num2str(trackinginfo(4)));
        set(handles.edit_lasttp,'String',num2str(trackinginfo(5)));
        set(handles.edit_templatecurrentframe,'String',num2str(trackinginfo(4)));
        set(handles.edit_signalcurrentframe,'String',num2str(trackinginfo(4)));
    case 5
        filetype = varargin{1};
        trackinginfo = varargin{2};
        imagelocation = varargin{3};
        channelnames = varargin{4};
        fileformat = varargin{5};%
        set(handles.edit_template,'String',num2str(trackinginfo(1)));
        set(handles.edit_nomin,'String',num2str(trackinginfo(2)));
        set(handles.edit_denomin,'String',num2str(trackinginfo(3)));
        set(handles.edit_firsttp,'String',num2str(trackinginfo(4)));
        set(handles.edit_lasttp,'String',num2str(trackinginfo(5)));
        set(handles.edit_templatecurrentframe,'String',num2str(trackinginfo(4)));
        set(handles.edit_signalcurrentframe,'String',num2str(trackinginfo(4)));
    case 12
        filetype = varargin{1};
        trackinginfo = varargin{2};
        imagelocation = varargin{3};
        channelnames = varargin{4};
        fileformat = varargin{5};
        currentcell = varargin{6};
        stageLoc = varargin{7};
        ind_cellpath = varargin{8};
        nucmask = varargin{9};
        cytomask = varargin{10};
        cellmask = varargin{11};
        selected_cells = varargin{12};
        set(handles.edit_template,'String',num2str(trackinginfo(1)));
        set(handles.edit_nomin,'String',num2str(trackinginfo(2)));
        set(handles.edit_denomin,'String',num2str(trackinginfo(3)));
        set(handles.edit_firsttp,'String',num2str(trackinginfo(4)));
        set(handles.edit_lasttp,'String',num2str(trackinginfo(5)));
        set(handles.edit_templatecurrentframe,'String',num2str(trackinginfo(4)));
        set(handles.edit_signalcurrentframe,'String',num2str(trackinginfo(4)));
        set(handles.edit_cellno,'String',num2str(currentcell));
        
        if ~isempty(ind_cellpath)
            handles.ind_cellpath = ind_cellpath;
        end
        if ~isempty(nucmask)  
            handles.currentnucmask = nucmask;
            clear nucmask;
        end
        if ~isempty(cytomask)  
            handles.currentcytomask = cytomask;
            clear cytomask;
        end
        if ~isempty(cellmask) 
            handles.currentcellmask = cellmask;
            clear cellmask;
        end
        
        if ~isempty(selected_cells)
            handles.selected_cells = selected_cells;
        end
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
        set(handles.radiobutton_petiffs,'Value',1);           
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
        [pathstr, ~, ~] = fileparts(fileformat);
        if isempty(pathstr)
            fileformat = fullfile(pwd,fileformat);
        else
            fileformat = fullfile(pathstr,fileformat);
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
        set(handles.radiobutton_customim,'Value',1);
        
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
         
        set(handles.edit_fileformat,'String',fileformat);
        handles.channelnames = channelnames;
        set(handles.edit_row,'String',row);
        set(handles.edit_col,'String',col);
        set(handles.edit_field,'String',field);
        set(handles.edit_plane,'String',plane);
        
        if strcmp(fileformat((end-2):end),'.nd')
            set(handles.edit_ndfilename,'String',fileformat);
            handles.ndfilename = fileformat;
            [notp stagePos stageName channelnames] = readndfile(fileformat);
            if notp==-1
                handles.initialframe = [];
                display([fileformat ' does not exist. Please re-define input ND file.']);
                handles.filetype = filetype;
                guidata(hObject, handles);
                return;
            end
            
            set(handles.edit_templatecurrentframe,'String',num2str(1));
            set(handles.edit_signalcurrentframe,'String',num2str(1));
            set(handles.popupmenu_stagePos,'String',stagePos);
            set(handles.popupmenu_stagePos,'Value',stageLoc);
            set(handles.edit_stageInfo,'String',stageName{stageLoc});
            handles.stageName = stageName;
            handles.channelnames = channelnames;
            prefix = fileformat(1:(end-3));
            handles.prefix = prefix;
            fileformat = [prefix '_%s_s' num2str(stageLoc) '_t%g.tif'];
            set(handles.edit_fileformat,'String',fileformat);
        end
        
end

handles.filetype = filetype;
set(handles.radiobutton_method1,'Value',1);
set(handles.togglebutton_nucleus,'Value',1);
set(handles.togglebutton_mean,'Value',1);
handles.cellrecogmethod = 1;
handles.segmenttype = 1;
handles.stattype = 1;
handles.greenflag=1;

guidata(hObject, handles);
handles = guidata(hObject);

minF = str2num(get(handles.edit_firsttp,'String'));
maxF = str2num(get(handles.edit_lasttp,'String'));

set(handles.slider_template,'Max',maxF);
set(handles.slider_template,'Min',minF);
set(handles.slider_template,'Value',str2num(get(handles.edit_templatecurrentframe,'String')));
set(handles.slider_template,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);

set(handles.slider_signal,'Max',maxF);
set(handles.slider_signal,'Min',minF);
set(handles.slider_signal,'Value',str2num(get(handles.edit_signalcurrentframe,'String')));
set(handles.slider_signal,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);
set(handles.edit_sourceframe,'String',num2str(minF) );
set(handles.edit_targetframes,'String',[num2str(minF) ':' num2str(maxF)]);


% save and update parameters
guidata(hObject, handles);

if length(varargin)==12
    
    guidata(hObject, handles);
    handles = guidata(hObject);
    pushbutton_load_Callback(hObject, eventdata, handles);
end




% UIWAIT makes MaskGenerating wait for user response (see UIRESUME)
% uiwait(handles.MaskGenerator);


% --- Outputs from this function are returned to the command line.
function varargout = MaskGenerating_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes during object creation, after setting all properties.
function slider_template_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_templatecurrentframe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_templatecurrentframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_templatecurrentframe as text
%        str2double(get(hObject,'String')) returns contents of edit_templatecurrentframe as a double
templatetp = str2num(get(hObject,'String'));
updatetemplate(handles);
set(handles.slider_template,'Value',templatetp);
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_templatecurrentframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_templatecurrentframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function slider_signal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function outputim = loadimage(filetype,fileformat,imlocation,tp,channelnames,pixelsreg,displayGate)
row = imlocation(1);
col = imlocation(2);
field = imlocation(3);
plane = imlocation(4);
channel1 = imlocation(5);
channel2 = imlocation(6);
totalCH = length(channelnames);

if channel2 == -1
    switch filetype
        case 1
            filename = sprintf(fileformat,row,col,field,plane,channel1,tp);
            outputim = imread(filename,'PixelRegion',pixelsreg);
        case 2
            outputim = imread(fileformat,'Index',totalCH*(tp-1)+channel1,'PixelRegion',pixelsreg);
        case 3
            filename = sprintf(fileformat,channelnames{channel1},tp);
            outputim = imread(filename,'PixelRegion',pixelsreg);
    end
    if isempty(displayGate)
        outputim = mat2gray(im2double(outputim));
    else
        ThresL = displayGate(1);
        ThresH = displayGate(2);
        outputim = mat2gray(im2double(outputim),[ThresL ThresH]);
    end
else
    switch filetype
        case 1
            filename = sprintf(fileformat,row,col,field,plane,channel1,tp);
            nominim = imread(filename,'PixelRegion',pixelsreg);
            filename = sprintf(fileformat,row,col,field,plane,channel2,tp);
            denomin = imread(filename,'PixelRegion',pixelsreg);
        case 2
            nominim = imread(fileformat,'Index',totalCH*(tp-1)+channel1,'PixelRegion',pixelsreg);
            denomin = imread(fileformat,'Index',totalCH*(tp-1)+channel2,'PixelRegion',pixelsreg);
        case 3
            filename = sprintf(fileformat,channelnames{channel1},tp);
            nominim = imread(filename,'PixelRegion',pixelsreg);
            filename = sprintf(fileformat,channelnames{channel2},tp);
            denomin = imread(filename,'PixelRegion',pixelsreg);
    end
    %ratioIm = calculateFRET(nominim,denomin);
    %outputim = mat2gray(ratioIm);
    outputim = mat2gray(im2double(nominim+2^9)./im2double(denomin+2^9),[1 1.5]);
end


function ratioIm = calculateFRET(nominIM,denomIM)
filterParam1 = 2;
filterParam2 = 2;

signalShiftN = 2^9;
signalShiftD = 2^9;

normN = double(ifft2(ifftshift(fftshift(fft2(double(im2int16(nominIM)))).*hbutter(nominIM,filterParam1,filterParam2))));
normD = double(ifft2(ifftshift(fftshift(fft2(double(im2int16(denomIM)))).*hbutter(denomIM,filterParam1,filterParam2))));

normN = normN+signalShiftN;
normD = normD+signalShiftD;

ratioIm =  normN./normD;


function iminfo = loadimageinfo(filetype,fileformat,imlocation,tp,channelnames,pixelsreg)
row = imlocation(1);
col = imlocation(2);
field = imlocation(3);
plane = imlocation(4);
channel = imlocation(5);
totalCH = length(channelnames);

switch filetype
    case 1
        filename = sprintf(fileformat,row,col,field,plane,channel,tp);
        iminfo = imfinfo(filename);
    case 2
        iminfo = imfinfo(fileformat,'Index',totalCH*(tp-1)+channel);
    case 3
        filename = sprintf(fileformat,channelnames{channel},tp);
        iminfo = imfinfo(filename);
end

function edit_signalcurrentframe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_signalcurrentframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_signalcurrentframe as text
%        str2double(get(hObject,'String')) returns contents of edit_signalcurrentframe as a double
signaltp = str2num(get(hObject,'String'));
updatesignal(handles);
set(handles.slider_signal,'Value',signaltp);
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_signalcurrentframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_signalcurrentframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
templateCH= str2num(get(handles.edit_template,'String'));
nominCH = str2num(get(handles.edit_nomin,'String'));
denomCH = str2num(get(handles.edit_denomin,'String'));
cellNo = str2num(get(handles.edit_cellno,'String'));
firsttp = str2num(get(handles.edit_firsttp,'String'));
lasttp = str2num(get(handles.edit_lasttp,'String'));

celltrackOUTfilename = ['celltrackOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(templateCH)];
maskOUTfilename = ['maskOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(templateCH)];


if exist([celltrackOUTfilename '.mat'],'file')
    load(celltrackOUTfilename);
else
    display([celltrackOUTfilename '.mat does not exist']);
    return;
end


if isempty(handles.currentnucmask)
    if exist([maskOUTfilename '.mat'],'file')
        
        nucName = ['nucmask_cell' num2str(cellNo)];
        cytoName = ['cytomask_cell' num2str(cellNo)];
        cellName = ['cellmask_cell' num2str(cellNo)];
        
        load(maskOUTfilename,nucName);
        load(maskOUTfilename,cytoName);
        load(maskOUTfilename,cellName);
        
        eval(['nucmask=' nucName ';']);
        eval(['cytomask=' cytoName ';']);
        eval(['cellmask=' cellName ';']);
        
        clear(nucName);
        clear(cytoName);
        clear(cellName);
    else
        display([maskOUTfilename '.mat does not exist']);
        return;
    end
    
else
    nucmask=handles.currentnucmask;
    cytomask=handles.currentcytomask;
    cellmask=handles.currentcellmask;
end

if exist('nucmask','var')
    cellsize = (size(nucmask(:,:,firsttp),1)-1)/2;
    set(handles.edit_cellsize,'String',num2str(cellsize));
    display(['changed cellsize to ' num2str(cellsize)]);
else
    cellsize = str2num(get(handles.edit_cellsize,'String'));
end


ThresL = str2num(get(handles.edit_gateL,'String'));
ThresH = str2num(get(handles.edit_gateH,'String'));
templateinfo = loadimageinfo(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane templateCH],firsttp,handles.channelnames);
imwidth = templateinfo.Width;
handles.imwidth = imwidth;
imheight = templateinfo.Height;
handles.imheight = imheight;
template = zeros(2*cellsize+1,2*cellsize+1,lasttp);
mysignal = zeros(2*cellsize+1,2*cellsize+1,lasttp);


ind_cellpath = handles.ind_cellpath; 


if ~isempty(handles.currentnucmask)
    if ~isempty( find( cellNo == handles.selected_cells ,1) )
        set(handles.togglebutton_selected,'Value',1);
    else
        set(handles.togglebutton_selected,'Value',0);
    end
    
    for tp = firsttp:lasttp
        if isempty(handles.ind_cellpath)
            ind_cellpath(tp,:) = cellpath{tp}(cellNo,:);
        end
        PosInd = find(cellpath{tp}(:,1)>0 & cellpath{tp}(:,2)>0)';
        ind_bg(tp,:) = round([mean(cellpath{tp}(PosInd,1)) mean(cellpath{tp}(PosInd,2))]);
        
        xL=max(ind_cellpath(tp,1)-cellsize,1);
        xR=min(ind_cellpath(tp,1)+cellsize,imwidth);
        yL=max(ind_cellpath(tp,2)-cellsize,1);
        yR=min(ind_cellpath(tp,2)+cellsize,imheight);
        
        
        if xR-xL == cellsize*2
            borderX = 1:(cellsize*2+1);
        elseif xR == imwidth
            shiftX = imwidth-xL;
            borderX = 1:(shiftX+1);
        elseif xL == 1
            shiftX = cellsize*2+1-xR;
            borderX = (xL+shiftX):(cellsize*2+1);
        end
        
        if yR-yL == cellsize*2
            borderY = 1:(cellsize*2+1);
        elseif yR == imheight
            shiftY = imheight-yL;
            borderY = 1:(shiftY+1);
        elseif yL == 1
            shiftY = cellsize*2+1-yR;
            borderY = (yL+shiftY):(cellsize*2+1);
        end
        
        
        template(borderY,borderX,tp) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane templateCH -1],tp,handles.channelnames,{[yL yR], [xL xR]},[]);
        mysignal(borderY,borderX,tp) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane nominCH denomCH],tp,handles.channelnames,{[yL yR], [xL xR]},[ThresL ThresH]);
        corners{tp} = [xL xR yL yR];
    end
    if ~isempty(handles.selected_cells)
        if ~isempty(find( cellNo == handles.selected_cells ))
            set(handles.togglebutton_selected,'Value',1);
        else
            set(handles.togglebutton_selected,'Value',0);
        end
    else
        set(handles.togglebutton_selected,'Value',0);
    end
else
    newnucmask = zeros(size(template,1),size(template,2),lasttp);
    newcellmask = zeros(size(template,1),size(template,2),lasttp);
    newcytomask = zeros(size(template,1),size(template,2),lasttp);
    if exist([maskOUTfilename '.mat'],'file')
        load(maskOUTfilename);
        clc;display(['Loaded maskes from ' maskOUTfilename]);
    end
    if isempty(handles.selected_cells)
        if exist('selected_cells','var')
            handles.selected_cells = selected_cells;
            if ~isempty(find( cellNo == selected_cells ))
                set(handles.togglebutton_selected,'Value',1);
            else
                set(handles.togglebutton_selected,'Value',0);
            end
        end
    end
    for tp = 1:size(ind_cellpath,1)
        
        if isempty(handles.ind_cellpath)
            ind_cellpath(tp,:) = cellpath{tp}(cellNo,:);
        end
        PosInd = find(cellpath{tp}(:,1)>0 & cellpath{tp}(:,2)>0)';
        ind_bg(tp,:) = round([mean(cellpath{tp}(PosInd,1)) mean(cellpath{tp}(PosInd,2))]);
        
        xL=max(ind_cellpath(tp,1)-cellsize,1);
        xR=min(ind_cellpath(tp,1)+cellsize,imwidth);
        yL=max(ind_cellpath(tp,2)-cellsize,1);
        yR=min(ind_cellpath(tp,2)+cellsize,imheight);
        
        if xR-xL == cellsize*2
            borderX = 1:(cellsize*2+1);
        elseif xR == imwidth
            shiftX = imwidth-xL;
            borderX = 1:(shiftX+1);
        elseif xL == 1
            shiftX = cellsize*2+1-xR;
            borderX = (xL+shiftX):(cellsize*2+1);
        end
        
        if yR-yL == cellsize*2
            borderY = 1:(cellsize*2+1);
        elseif yR == imheight
            shiftY = imheight-yL;
            borderY = 1:(shiftY+1);
        elseif yL == 1
            shiftY = cellsize*2+1-yR;
            borderY = (yL+shiftY):(cellsize*2+1);
        end
        
        template(borderY,borderX,tp) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane templateCH -1],tp,handles.channelnames,{[yL yR], [xL xR]},[]);
        mysignal(borderY,borderX,tp) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane nominCH denomCH],tp,handles.channelnames,{[yL yR], [xL xR]},[ThresL ThresH]);
        corners{tp} = [xL xR yL yR];
        switch handles.cellrecogmethod
            case 1
                [BW] = templateToCentroid(template(:,:,tp),handles,tp);
            case 2
                [BW] = templateToCentroid2(template(:,:,tp),handles,tp);
        end
        if exist('nucmask','var') && size(nucmask,3)>=tp  && ~isempty(find(nucmask(:,:,tp)))
            newnucmask(:,:,tp) = nucmask(:,:,tp);
        else
            newnucmask(:,:,tp) = BW;
        end
        
        if exist('cytomask','var') && size(cytomask,3)>=tp  && ~isempty(find(cytomask(:,:,tp)))
            newcytomask(:,:,tp) = cytomask(:,:,tp);
        end
        
        if exist('cellmask','var')  && size(cellmask,3)>=tp  && ~isempty(find(cellmask(:,:,tp)))
            newcellmask(:,:,tp) = cellmask(:,:,tp);
        end
        
        
    end
    handles.currentnucmask = newnucmask;
    handles.currentcytomask = newcytomask;
    handles.currentcellmask = newcellmask;
    
end

handles.ind_cellpath = ind_cellpath;
handles.ind_bg = ind_bg;
handles.template = template;
handles.mysignal = mysignal;
handles.corners = corners;
guidata(hObject, handles);
handles = guidata(hObject);

updatexyplots(handles);
updatesignalplot(handles);
updatetemplate(handles);
updatesignal(handles);
minF = str2num(get(handles.edit_firsttp,'String'));
maxF = str2num(get(handles.edit_lasttp,'String'));

set(handles.slider_template,'Max',maxF);
set(handles.slider_template,'Min',minF);
set(handles.slider_template,'Value',str2num(get(handles.edit_templatecurrentframe,'String')));
set(handles.slider_template,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);

set(handles.slider_signal,'Max',maxF);
set(handles.slider_signal,'Min',minF);
set(handles.slider_signal,'Value',str2num(get(handles.edit_signalcurrentframe,'String')));
set(handles.slider_signal,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);

guidata(hObject, handles);



function updatesignalplot(handles)
nucmask  = handles.currentnucmask;
cellmask = handles.currentcellmask;
cytomask = handles.currentcytomask;
ThresL = str2num(get(handles.edit_gateL,'String'));
ThresH = str2num(get(handles.edit_gateH,'String'));

signalIm = (handles.mysignal*(ThresH-ThresL))+ThresL;
firsttp = str2num(get(handles.edit_firsttp,'String'));
lasttp = str2num(get(handles.edit_lasttp,'String'));

sindex=1;
plotsignal = [];
for tp = firsttp:lasttp
    
    switch handles.segmenttype
        case 1
            BW = nucmask(:,:,tp);
            [Rs Cs] = find(BW);
            if ~isempty(find(BW))
                finalSignal = signalIm(Rs,Cs,tp);
                
                switch handles.stattype
                    case 1
                        plotsignal(sindex,:) = [tp mean(finalSignal(:))];
                    case 2
                        plotsignal(sindex,:) = [tp median(finalSignal(:))];
                    case 3
                        plotsignal(sindex,:) = [tp sum(finalSignal(:))];
                end
                sindex = sindex+1;
            end
            
        case 2
            if ~isempty(find(cytomask(:,:,tp)))
                BW = cytomask(:,:,tp);
            else
                BW = cellmask(:,:,tp) & ~ nucmask(:,:,tp);
            end
            [Rs Cs] = find(BW);
            if ~isempty(find(BW))
                finalSignal = signalIm(Rs,Cs,tp);
                
                switch handles.stattype
                    case 1
                        plotsignal(sindex,:) = [tp mean(finalSignal(:))];
                    case 2
                        plotsignal(sindex,:) = [tp median(finalSignal(:))];
                    case 3
                        plotsignal(sindex,:) = [tp sum(finalSignal(:))];
                end
                sindex = sindex+1;
            end
        case 3
            BW1 = nucmask(:,:,tp);
            
            if ~isempty(find(cytomask(:,:,tp)))
                BW2 = cytomask(:,:,tp);
            else
                BW2 = cellmask(:,:,tp) & ~ nucmask(:,:,tp);
            end
            
            [Rs1 Cs1] = find(BW1);
            [Rs2 Cs2] = find(BW2);
            if ~isempty(find(BW1)) & ~isempty(find(BW2))
                nuclear = signalIm(Rs1,Cs1,tp);
                cytosol = signalIm(Rs2,Cs2,tp);
                
                switch handles.stattype
                    case 1
                        plotsignal(sindex,:) = [tp mean(nuclear(:))/mean(cytosol(:))];
                    case 2
                        plotsignal(sindex,:) = [tp median(nuclear(:))/median(cytosol(:))];
                    case 3
                        plotsignal(sindex,:) = [tp sum(nuclear(:))/sum(cytosol(:))];
                end
                
                sindex = sindex+1;
            end
    end
end

axes(handles.axes_SignalTime);
if ~isempty(plotsignal)
    plot(plotsignal(:,1),plotsignal(:,2),'b.');
end
minF = str2num(get(handles.edit_firsttp,'String'));
maxF = str2num(get(handles.edit_lasttp,'String'));
xlim([minF maxF]);


% --- Executes on slider movement.
function slider_template_Callback(hObject, eventdata, handles)
% hObject    handle to slider_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
templatetp = round(get(hObject,'Value'));
set(handles.edit_templatecurrentframe,'String',num2str(templatetp));
guidata(hObject, handles);
updatetemplate(handles);

% --- Executes on slider movement.
function slider_signal_Callback(hObject, eventdata, handles)
% hObject    handle to slider_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
signaltp = round(get(hObject,'Value'));
set(handles.edit_signalcurrentframe,'String',num2str(signaltp));
guidata(hObject, handles);
updatesignal(handles);

function updatexyplots(handles)
ind_cellpath = handles.ind_cellpath;
ind_bg = handles.ind_bg;
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
templateCH= str2num(get(handles.edit_template,'String'));
nominCH = str2num(get(handles.edit_nomin,'String'));
denomCH = str2num(get(handles.edit_denomin,'String'));
cellNo = str2num(get(handles.edit_cellno,'String'));
cellsize = str2num(get(handles.edit_cellsize,'String'));
firsttp = str2num(get(handles.edit_firsttp,'String'));
lasttp = str2num(get(handles.edit_lasttp,'String'));

axes(handles.axes_XY);
plot(ind_cellpath(firsttp:lasttp,1)-ind_bg(firsttp:lasttp,1),ind_cellpath(firsttp:lasttp,2)-ind_bg(firsttp:lasttp,2),'b.-');hold on;
plot(ind_cellpath(firsttp,1)-ind_bg(firsttp,1),ind_cellpath(firsttp,2)-ind_bg(firsttp,2),'og','MarkerFaceColor','g');
plot(ind_cellpath(lasttp,1)-ind_bg(lasttp,1),ind_cellpath(lasttp,2)-ind_bg(lasttp,2),'or','MarkerFaceColor','r');
xlabel('X');
ylabel('Y');
hold off;

axes(handles.axes_XTime);
plot(firsttp:lasttp,ind_cellpath(firsttp:lasttp,1)-ind_bg(firsttp:lasttp,1),'b.-');hold on;
plot(firsttp,ind_cellpath(firsttp,1)-ind_bg(firsttp,1),'og','MarkerFaceColor','g');
plot(lasttp,ind_cellpath(lasttp,1)-ind_bg(lasttp,1),'or','MarkerFaceColor','r');
xlim([firsttp lasttp]);
hold off;

axes(handles.axes_YTime); 
plot(firsttp:lasttp,ind_cellpath(firsttp:lasttp,2)-ind_bg(firsttp:lasttp,2),'b.-');hold on;
plot(firsttp,ind_cellpath(firsttp,2)-ind_bg(firsttp,2),'og','MarkerFaceColor','g');
plot(lasttp,ind_cellpath(lasttp,2)-ind_bg(lasttp,2),'or','MarkerFaceColor','r');
xlim([firsttp lasttp]);
hold off;


% --- Executes on button press in pushbutton_roisubstract.
function pushbutton_roisubstract_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_roisubstract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tp = str2num(get(handles.edit_templatecurrentframe,'String'));
BW = handles.currentnucmask(:,:,tp);
h = imfreehand(handles.axes_template);
accepted_pos = wait(h);
sBW = createMask(h);
delete(h);
finalBW = BW&~sBW;
handles.currentnucmask(:,:,tp) = finalBW;
guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);

% --- Executes on button press in pushbutton_roiadd.
function pushbutton_roiadd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_roiadd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tp = str2num(get(handles.edit_templatecurrentframe,'String'));
BW = handles.currentnucmask(:,:,tp);
h = imfreehand(handles.axes_template);
accepted_pos = wait(h);
sBW = createMask(h);
delete(h);
finalBW = BW|sBW;
handles.currentnucmask(:,:,tp) = finalBW;
guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);


% --- Executes on button press in pushbutton_centerxy.
function pushbutton_centerxy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_centerxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
templateCH= str2num(get(handles.edit_template,'String'));
nominCH = str2num(get(handles.edit_nomin,'String'));
denomCH = str2num(get(handles.edit_denomin,'String'));
cellNo = str2num(get(handles.edit_cellno,'String'));
ThresL = str2num(get(handles.edit_gateL,'String'));
ThresH = str2num(get(handles.edit_gateH,'String'));
ind_cellpath = handles.ind_cellpath;
ind_bg = handles.ind_bg;
cellsize = str2num(get(handles.edit_cellsize,'String'));
tp = str2num(get(handles.edit_templatecurrentframe,'String'));

BW = handles.currentnucmask(:,:,tp);
S  = regionprops(BW, 'centroid');

if isempty(find(BW==0)) | isempty(find(BW==1))
    display('Nuclear mask not compatible');
    return;

else
    x = round(S.Centroid(1));
    y = round(S.Centroid(2));
end

xg = cellsize+1;
yg = cellsize+1;
if abs(yg-y)>1 | abs(xg-x)>1
    oldcorners = handles.corners{tp};
    ind_cellpath(tp,:) = [x+oldcorners(1) y+oldcorners(3)];
    
    xL=max(ind_cellpath(tp,1)-cellsize,1);
    xR=min(ind_cellpath(tp,1)+cellsize,handles.imwidth);
    yL=max(ind_cellpath(tp,2)-cellsize,1);
    yR=min(ind_cellpath(tp,2)+cellsize,handles.imheight);
    
    if xR-xL == cellsize*2
        borderX = 1:(cellsize*2+1);
    elseif xR == handles.imwidth
        shiftX = handles.imwidth-xL;
        borderX = 1:(shiftX+1);
    elseif xL == 1
        shiftX = cellsize*2+1-xR;
        borderX = (xL+shiftX):(cellsize*2+1);
    end
    
    if yR-yL == cellsize*2
        borderY = 1:(cellsize*2+1);
    elseif yR == handles.imheight
        shiftY = handles.imheight-yL;
        borderY = 1:(shiftY+1);
    elseif yL == 1
        shiftY = cellsize*2+1-yR;
        borderY = (yL+shiftY):(cellsize*2+1);
    end
    
    template(borderY,borderX) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane templateCH -1],tp,handles.channelnames,{[yL yR], [xL xR]},[]);
    signalIm(borderY,borderX) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane nominCH denomCH],tp,handles.channelnames,{[yL yR], [xL xR]},[ThresL ThresH]);
    
    handles.ind_cellpath = ind_cellpath;
    handles.template(:,:,tp) = template;
    handles.mysignal(:,:,tp) = signalIm;
    handles.corners{tp} = [xL xR yL yR];
    
    [oldR oldC] = find(BW);
    
    newR = round(oldR+yg-y);
    newC = round(oldC+xg-x);
    newBW = zeros(size(BW,1),size(BW,2));
    for i=1:length(newR)
        newBW(newR(i),newC(i)) = 1;
    end

    BW2 = handles.currentcellmask(:,:,tp);
    if isempty(find(BW2))
        newBW2 = zeros(size(BW2,1),size(BW2,2));
    else
        [oldR oldC] = find(BW2);
        newR = round(oldR+yg-y);
        newC = round(oldC+xg-x);
        newBW2 = zeros(size(BW2,1),size(BW2,2));
        for i=1:length(newR)
            newBW2(newR(i),newC(i)) = 1;
        end
        
    end
    
    BW3 = handles.currentcytomask(:,:,tp);
    if isempty(find(BW3))
        newBW3 = zeros(size(BW3,1),size(BW3,2));
    else
        [oldR oldC] = find(BW3);
        newR = round(oldR+yg-y);
        newC = round(oldC+xg-x);
        newBW3 = zeros(size(BW3,1),size(BW3,2));
        for i=1:length(newR)
            newBW3(newR(i),newC(i)) = 1;
        end
        
    end
    
    
    

else
    newBW = handles.currentnucmask(:,:,tp);
    newBW2 = handles.currentcellmask(:,:,tp);
    newBW3 = handles.currentcytomask(:,:,tp);
end

currentnucmask = handles.currentnucmask;
currentnucmask(:,:,tp) = newBW;
handles.currentnucmask = currentnucmask;

currentcellmask = handles.currentcellmask;
currentcellmask(:,:,tp) = newBW2;
handles.currentcellmask = currentcellmask;


currentcytomask = handles.currentcytomask;
currentcytomask(:,:,tp) = newBW3;
handles.currentcytomask = currentcytomask;

guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);
updatexyplots(handles);
updatesignalplot(handles);

function updatetemplate(handles)
tp = str2num(get(handles.edit_templatecurrentframe,'String'));
template = handles.template(:,:,tp);
nuc = handles.currentnucmask(:,:,tp);
cyto = handles.currentcytomask(:,:,tp);
cell = handles.currentcellmask(:,:,tp);

corners = handles.corners{tp};
xL = corners(1);
xR = corners(2);
yL = corners(3);
yR = corners(4);
ind_cellpath = handles.ind_cellpath;

bwNucEdge = bwperim(nuc);
bwCytoEdge = bwperim(cyto);
bwCellEdge = bwperim(cell);
cp = zeros(size(bwNucEdge));
cellsize = (size(nuc,1)-1)/2;
cp(cellsize,cellsize) = 1;

bwFinal = bwNucEdge| bwCellEdge | cp | bwCytoEdge;

imAdj = template-bwFinal;
imOut=cat(3,max(imAdj,(bwNucEdge | cp)),max(imAdj,bwCytoEdge),max(imAdj,bwCellEdge));
imshow(imOut,[],'Parent',handles.axes_template);

% set(handles.axes_template,'NextPlot','add');
% 
% plot(handles.axes_template,x,y,'.g','MarkerSize',1);
% plot(handles.axes_template,,'xr');
% set(handles.axes_template,'NextPlot','replace');

function updatesignal(handles)
tp = str2num(get(handles.edit_signalcurrentframe,'String'));
mysignal = handles.mysignal(:,:,tp);
nuc = handles.currentnucmask(:,:,tp);
cyto = handles.currentcytomask(:,:,tp);
cell = handles.currentcellmask(:,:,tp);

corners = handles.corners{tp};
xL = corners(1);
xR = corners(2);
yL = corners(3);
yR = corners(4);
ind_cellpath = handles.ind_cellpath;

bwNucEdge = bwperim(nuc);
bwCytoEdge = bwperim(cyto);
bwCellEdge = bwperim(cell);
cp = zeros(size(bwNucEdge));
cellsize = (size(nuc,1)-1)/2;
cp(cellsize,cellsize) = 1;

bwFinal = bwNucEdge| bwCellEdge | cp | bwCytoEdge;
imAdj = mysignal-bwFinal;
imOut=cat(3,max(imAdj,(bwNucEdge | cp)),max(imAdj,bwCytoEdge),max(imAdj,bwCellEdge));
imshow(imOut,[0 1],'Parent',handles.axes_signal);


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
templateCH= str2num(get(handles.edit_template,'String'));
cellNo = str2num(get(handles.edit_cellno,'String'));
firsttp = str2num(get(handles.edit_firsttp,'String'));
lasttp = str2num(get(handles.edit_lasttp,'String'));

celltrackOUTfilename = ['celltrackOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(templateCH)];
maskOUTfilename = ['maskOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(templateCH)];
load(celltrackOUTfilename);
if exist([maskOUTfilename '.mat'],'file')
    load(maskOUTfilename);
end

for tp = firsttp:lasttp
    cellpath{tp}(cellNo,:) = handles.ind_cellpath(tp,:);
end
nucmask{cellNo}  = handles.currentnucmask;
cytomask{cellNo} = handles.currentcytomask;
cellmask{cellNo} = handles.currentcellmask;

selected_cells = handles.selected_cells;
save(celltrackOUTfilename,'cellpath','sisterList','bg')
save(maskOUTfilename,'nucmask','cellmask','cytomask','selected_cells'); 


% --- Executes on button press in pushbutton_adoptnuc.
function pushbutton_adoptnuc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_adoptnuc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sourceF = str2num(get(handles.edit_sourceframe,'String'));
targetF = str2num(get(handles.edit_targetframes,'String'));

nucmask  = handles.currentnucmask;

if length(sourceF)==1 & length(targetF)>0
    for i=targetF
        nucmask(:,:,i) = nucmask(:,:,sourceF);
    end

    handles.currentnucmask = nucmask;

    guidata(hObject, handles);
    handles = guidata(hObject);
    updatexyplots(handles);
    updatesignalplot(handles);
    updatetemplate(handles);
    updatesignal(handles);
else
    display('Size of source and/or target frames are not appropriate.');
    return;
end

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_pausetemplate.
function pushbutton_pausetemplate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pausetemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.greenflag = 0;
guidata(hObject, handles);

% --- Executes on button press in pushbutton_playtemplate.
function pushbutton_playtemplate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_playtemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currenttp = str2num(get(handles.edit_templatecurrentframe,'String'));
lasttp = str2num(get(handles.edit_lasttp,'String'));
handles.greenflag = 1;
guidata(hObject, handles);
handles = guidata(hObject);
for tp = currenttp:lasttp
    set(handles.edit_templatecurrentframe,'String',num2str(tp));
    set(handles.slider_template,'Value',tp);
    updatetemplate(handles);
    pause(1/str2num(get(handles.edit_fps,'String')));
    handles = guidata(hObject);
    if handles.greenflag==0
        break;
    end
end
handles.greenflag = 1;
guidata(hObject, handles);

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_lastframetemplate.
function pushbutton_lastframetemplate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_lastframetemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lasttp = str2num(get(handles.edit_lasttp,'String'));

set(handles.edit_templatecurrentframe,'String',num2str(lasttp));
set(handles.slider_template,'Value',lasttp);
updatetemplate(handles);

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_firstframettemplate.
function pushbutton_firstframettemplate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_firstframettemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

firsttp = str2num(get(handles.edit_firsttp,'String'));

set(handles.edit_templatecurrentframe,'String',num2str(firsttp));
set(handles.slider_template,'Value',firsttp);
updatetemplate(handles);

% --- Executes on button press in togglebutton_nucleus.
function togglebutton_nucleus_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_nucleus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton_cytosol.
function togglebutton_cytosol_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_cytosol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton_nucPcyto.
function togglebutton_nucPcyto_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_nucPcyto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton_integrated.
function togglebutton_integrated_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_integrated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton_median.
function togglebutton_median_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton_mean.
function togglebutton_mean_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_invertedLog.
function checkbox_invertedLog_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_invertedLog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_invertedLog



function edit_displayHThres_Callback(hObject, eventdata, handles)
% hObject    handle to edit_displayHThres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_displayHThres as text
%        str2double(get(hObject,'String')) returns contents of edit_displayHThres as a double


% --- Executes during object creation, after setting all properties.
function edit_displayHThres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_displayHThres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_displayLThres_Callback(hObject, eventdata, handles)
% hObject    handle to edit_displayLThres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_displayLThres as text
%        str2double(get(hObject,'String')) returns contents of edit_displayLThres as a double


% --- Executes during object creation, after setting all properties.
function edit_displayLThres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_displayLThres (see GCBO)
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



function edit_template_Callback(hObject, eventdata, handles)
% hObject    handle to edit_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_template as text
%        str2double(get(hObject,'String')) returns contents of edit_template as a double


% --- Executes during object creation, after setting all properties.
function edit_template_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_template (see GCBO)
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



function edit_cellno_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cellno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cellno as text
%        str2double(get(hObject,'String')) returns contents of edit_cellno as a double


% --- Executes during object creation, after setting all properties.
function edit_cellno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cellno (see GCBO)
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


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_nomin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nomin as text
%        str2double(get(hObject,'String')) returns contents of edit_nomin as a double


% --- Executes during object creation, after setting all properties.
function edit_nomin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_denomin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_denomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_denomin as text
%        str2double(get(hObject,'String')) returns contents of edit_denomin as a double


% --- Executes during object creation, after setting all properties.
function edit_denomin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_denomin (see GCBO)
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


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_firsttp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_firsttp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_firsttp as text
%        str2double(get(hObject,'String')) returns contents of edit_firsttp as a double

minF = str2num(get(hObject,'String'));
maxF = str2num(get(handles.edit_lasttp,'String'));

set(handles.slider_template,'Max',maxF);
set(handles.slider_template,'Min',minF);
set(handles.slider_template,'Value',minF);
set(handles.slider_template,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);

set(handles.slider_signal,'Max',maxF);
set(handles.slider_signal,'Min',minF);
set(handles.slider_signal,'Value',minF);
set(handles.slider_signal,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_firsttp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_firsttp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_lasttp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lasttp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lasttp as text
%        str2double(get(hObject,'String')) returns contents of edit_lasttp as a double

minF = str2num(get(handles.edit_firsttp,'String'));
maxF = str2num(get(hObject,'String'));

set(handles.slider_template,'Max',maxF);
set(handles.slider_template,'Min',minF);
set(handles.slider_template,'Value',str2num(get(handles.edit_templatecurrentframe,'String')));
set(handles.slider_template,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);

set(handles.slider_signal,'Max',maxF);
set(handles.slider_signal,'Min',minF);
set(handles.slider_signal,'Value',str2num(get(handles.edit_signalcurrentframe,'String')));
set(handles.slider_signal,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_lasttp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lasttp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton_definenuc.
function pushbutton_definenuc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_definenuc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tp = str2num(get(handles.edit_templatecurrentframe,'String'));
h = imfreehand(handles.axes_template);
accepted_pos = wait(h);
BW = createMask(h);
delete(h);
handles.currentnucmask(:,:,tp) = BW;
guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);
updatesignalplot(handles);

% --- Executes on button press in pushbutton_definecell.
function pushbutton_definecell_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_definecell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tp = str2num(get(handles.edit_templatecurrentframe,'String'));
h = imfreehand(handles.axes_template);
accepted_pos = wait(h);
BW = createMask(h);
delete(h);
handles.currentcellmask(:,:,tp) = BW;
guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);
updatesignalplot(handles);

% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_sourceframe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sourceframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sourceframe as text
%        str2double(get(hObject,'String')) returns contents of edit_sourceframe as a double


% --- Executes during object creation, after setting all properties.
function edit_sourceframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sourceframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_targetframes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_targetframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_targetframes as text
%        str2double(get(hObject,'String')) returns contents of edit_targetframes as a double


% --- Executes during object creation, after setting all properties.
function edit_targetframes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_targetframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_storenuc.
function pushbutton_storenuc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_storenuc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)






function BW = templateToCentroid(M,handles,tp)
tp = str2num(get(handles.edit_templatecurrentframe,'String'));
xg = round(size(M,2)/2);
yg = round(size(M,1)/2);
maxI = str2num(get(handles.edit_maxI,'String'));
invertedLog = get(handles.checkbox_invertedLog,'Value');

if invertedLog
    inverted = ((maxI)-1)-M;
end
BWc = zeros(size(M));
for i=2:0.7:6
    edgedIm = edge(M,'canny',0,i);
    BW = imfill(edgedIm,'holes');
    
    BW = bwmorph(BW,'open',1);
    BW = bwselect(BW,xg,yg);
    
    BWc = BWc | BW;

end
BW=BWc;

function outBW = templateToCentroid2(M,handles,tp)

M = handles.template(:,:,tp);
xg = round(size(M,2)/2);
yg = round(size(M,1)/2);
maxI = str2num(get(handles.edit_maxI,'String'));
invertedLog = get(handles.checkbox_invertedLog,'Value');
outerbox = str2num(get(handles.edit_cellsize,'String'));

if invertedLog
    M = ((maxI)-1)-M;
end

M(1,:) = 0;
M(end,:) = 0;
M(:,1) = 0;
M(:,end) = 0;
outBW = modchenvese(M,150,0.1,outerbox);
outBW(1,:) = 0;
outBW(end,:) = 0;
outBW(:,1) = 0;
outBW(:,end) = 0;


%   Adapted from code by Yue Wu (yue.wu@tufts.edu)
%   http://www.mathworks.com/matlabcentral/fileexchange/23445
function seg = modchenvese(I,num_iter,mu,outerbox)
I = imadjust(I);
s = outerbox./min(size(I,1),size(I,2)); % resize scale

n = zeros(size(I));
midY = round(size(n,1)/2);
midX = round(size(n,2)/2);
boxsize = round(outerbox/2*.9);
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
    phi_{j} = phi(:,:,j);
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
    if Q<=dt*.18^2
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
H(idx1)=1;
for i=1:length(idx2)
    H(idx2(i))=1/2*(1+z(idx2(i))/Epsilon+1/pi*sin(pi*z(idx2(i))/Epsilon));
end;
    
    
    
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
        handles.cellrecogmethod = 1; 
    case 'radiobutton_method2'
        handles.cellrecogmethod = 2;
    otherwise
end
guidata(hObject, handles);  


% --- Executes on button press in pushbutton_storecyto.
function pushbutton_storecyto_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_storecyto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in uipanel_segment.
function uipanel_segment_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_segment 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'togglebutton_nucleus'
        handles.segmenttype = 1; 
    case 'togglebutton_cytosol'
        handles.segmenttype = 2;
    case 'togglebutton_nucPcyto'
        handles.segmenttype = 3;
end
guidata(hObject, handles);  
handles = guidata(hObject);  
updatesignalplot(handles);

% --- Executes when selected object is changed in uipanel_stat.
function uipanel_stat_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_stat 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'togglebutton_mean'
        handles.stattype = 1; 
    case 'togglebutton_median'
        handles.stattype = 2;
    case 'togglebutton_integrated'
        handles.stattype = 3;
end
guidata(hObject, handles);  
handles = guidata(hObject);  
updatesignalplot(handles);


% --- Executes on button press in pushbutton_shiftright.
function pushbutton_shiftright_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_shiftright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
templateCH= str2num(get(handles.edit_template,'String'));
nominCH = str2num(get(handles.edit_nomin,'String'));
denomCH = str2num(get(handles.edit_denomin,'String'));
cellNo = str2num(get(handles.edit_cellno,'String'));
ThresL = str2num(get(handles.edit_gateL,'String'));
ThresH = str2num(get(handles.edit_gateH,'String'));
ind_cellpath = handles.ind_cellpath;
ind_bg = handles.ind_bg;
cellsize = str2num(get(handles.edit_cellsize,'String'));
tp = str2num(get(handles.edit_templatecurrentframe,'String'));

oldcorners = handles.corners{tp};

ind_cellpath(tp,:) = [ind_cellpath(tp,1)+ str2num(get(handles.edit_shiftsize,'String')) ind_cellpath(tp,2)];

xL=max(ind_cellpath(tp,1)-cellsize,1);
xR=min(ind_cellpath(tp,1)+cellsize,handles.imwidth);
yL=max(ind_cellpath(tp,2)-cellsize,1);
yR=min(ind_cellpath(tp,2)+cellsize,handles.imheight);

if xR-xL == cellsize*2
    borderX = 1:(cellsize*2+1);
elseif xR == handles.imwidth
    shiftX = handles.imwidth-xL;
    borderX = 1:(shiftX+1);
elseif xL == 1
    shiftX = cellsize*2+1-xR;
    borderX = (xL+shiftX):(cellsize*2+1);
end

if yR-yL == cellsize*2
    borderY = 1:(cellsize*2+1);
elseif yR == handles.imheight
    shiftY = handles.imheight-yL;
    borderY = 1:(shiftY+1);
elseif yL == 1
    shiftY = cellsize*2+1-yR;
    borderY = (yL+shiftY):(cellsize*2+1);
end

template(borderY,borderX) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane templateCH -1],tp,handles.channelnames,{[yL yR], [xL xR]},[]);
signalIm(borderY,borderX) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane nominCH denomCH],tp,handles.channelnames,{[yL yR], [xL xR]},[ThresL ThresH]);

handles.ind_cellpath = ind_cellpath;
handles.template(:,:,tp) = template;
handles.mysignal(:,:,tp) = signalIm;
handles.corners{tp} = [xL xR yL yR];

guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);
updatexyplots(handles);
updatesignalplot(handles);

% --- Executes on button press in pushbutton_shiftup.
function pushbutton_shiftup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_shiftup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
templateCH= str2num(get(handles.edit_template,'String'));
nominCH = str2num(get(handles.edit_nomin,'String'));
denomCH = str2num(get(handles.edit_denomin,'String'));
cellNo = str2num(get(handles.edit_cellno,'String'));
ThresL = str2num(get(handles.edit_gateL,'String'));
ThresH = str2num(get(handles.edit_gateH,'String'));
ind_cellpath = handles.ind_cellpath;
ind_bg = handles.ind_bg;
cellsize = str2num(get(handles.edit_cellsize,'String'));
tp = str2num(get(handles.edit_templatecurrentframe,'String'));

oldcorners = handles.corners{tp};

ind_cellpath(tp,:) = [ind_cellpath(tp,1) ind_cellpath(tp,2)- str2num(get(handles.edit_shiftsize,'String'))];

xL=max(ind_cellpath(tp,1)-cellsize,1);
xR=min(ind_cellpath(tp,1)+cellsize,handles.imwidth);
yL=max(ind_cellpath(tp,2)-cellsize,1);
yR=min(ind_cellpath(tp,2)+cellsize,handles.imheight);

if xR-xL == cellsize*2
    borderX = 1:(cellsize*2+1);
elseif xR == handles.imwidth
    shiftX = handles.imwidth-xL;
    borderX = 1:(shiftX+1);
elseif xL == 1
    shiftX = cellsize*2+1-xR;
    borderX = (xL+shiftX):(cellsize*2+1);
end

if yR-yL == cellsize*2
    borderY = 1:(cellsize*2+1);
elseif yR == handles.imheight
    shiftY = handles.imheight-yL;
    borderY = 1:(shiftY+1);
elseif yL == 1
    shiftY = cellsize*2+1-yR;
    borderY = (yL+shiftY):(cellsize*2+1);
end

template(borderY,borderX) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane templateCH -1],tp,handles.channelnames,{[yL yR], [xL xR]},[]);
signalIm(borderY,borderX) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane nominCH denomCH],tp,handles.channelnames,{[yL yR], [xL xR]},[ThresL ThresH]);

handles.ind_cellpath = ind_cellpath;
handles.template(:,:,tp) = template;
handles.mysignal(:,:,tp) = signalIm;
handles.corners{tp} = [xL xR yL yR];

guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);
updatexyplots(handles);
updatesignalplot(handles);

% --- Executes on button press in pushbutton_shiftleft.
function pushbutton_shiftleft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_shiftleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
templateCH= str2num(get(handles.edit_template,'String'));
nominCH = str2num(get(handles.edit_nomin,'String'));
denomCH = str2num(get(handles.edit_denomin,'String'));
cellNo = str2num(get(handles.edit_cellno,'String'));
ThresL = str2num(get(handles.edit_gateL,'String'));
ThresH = str2num(get(handles.edit_gateH,'String'));
ind_cellpath = handles.ind_cellpath;
ind_bg = handles.ind_bg;
cellsize = str2num(get(handles.edit_cellsize,'String'));
tp = str2num(get(handles.edit_templatecurrentframe,'String'));

oldcorners = handles.corners{tp};

ind_cellpath(tp,:) = [ind_cellpath(tp,1)- str2num(get(handles.edit_shiftsize,'String')) ind_cellpath(tp,2)];

xL=max(ind_cellpath(tp,1)-cellsize,1);
xR=min(ind_cellpath(tp,1)+cellsize,handles.imwidth);
yL=max(ind_cellpath(tp,2)-cellsize,1);
yR=min(ind_cellpath(tp,2)+cellsize,handles.imheight);

if xR-xL == cellsize*2
    borderX = 1:(cellsize*2+1);
elseif xR == handles.imwidth
    shiftX = handles.imwidth-xL;
    borderX = 1:(shiftX+1);
elseif xL == 1
    shiftX = cellsize*2+1-xR;
    borderX = (xL+shiftX):(cellsize*2+1);
end

if yR-yL == cellsize*2
    borderY = 1:(cellsize*2+1);
elseif yR == handles.imheight
    shiftY = handles.imheight-yL;
    borderY = 1:(shiftY+1);
elseif yL == 1
    shiftY = cellsize*2+1-yR;
    borderY = (yL+shiftY):(cellsize*2+1);
end

template(borderY,borderX) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane templateCH -1],tp,handles.channelnames,{[yL yR], [xL xR]},[]);
signalIm(borderY,borderX) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane nominCH denomCH],tp,handles.channelnames,{[yL yR], [xL xR]},[ThresL ThresH]);

handles.ind_cellpath = ind_cellpath;
handles.template(:,:,tp) = template;
handles.mysignal(:,:,tp) = signalIm;
handles.corners{tp} = [xL xR yL yR];

guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);
updatexyplots(handles);
updatesignalplot(handles);

% --- Executes on button press in pushbutton_shiftdown.
function pushbutton_shiftdown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_shiftdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
templateCH= str2num(get(handles.edit_template,'String'));
nominCH = str2num(get(handles.edit_nomin,'String'));
denomCH = str2num(get(handles.edit_denomin,'String'));
cellNo = str2num(get(handles.edit_cellno,'String'));
ThresL = str2num(get(handles.edit_gateL,'String'));
ThresH = str2num(get(handles.edit_gateH,'String'));

ind_cellpath = handles.ind_cellpath;
ind_bg = handles.ind_bg;
cellsize = str2num(get(handles.edit_cellsize,'String'));
tp = str2num(get(handles.edit_templatecurrentframe,'String'));

oldcorners = handles.corners{tp};

ind_cellpath(tp,:) = [ind_cellpath(tp,1) ind_cellpath(tp,2)+ str2num(get(handles.edit_shiftsize,'String'))];

xL=max(ind_cellpath(tp,1)-cellsize,1);
xR=min(ind_cellpath(tp,1)+cellsize,handles.imwidth);
yL=max(ind_cellpath(tp,2)-cellsize,1);
yR=min(ind_cellpath(tp,2)+cellsize,handles.imheight);

if xR-xL == cellsize*2
    borderX = 1:(cellsize*2+1);
elseif xR == handles.imwidth
    shiftX = handles.imwidth-xL;
    borderX = 1:(shiftX+1);
elseif xL == 1
    shiftX = cellsize*2+1-xR;
    borderX = (xL+shiftX):(cellsize*2+1);
end

if yR-yL == cellsize*2
    borderY = 1:(cellsize*2+1);
elseif yR == handles.imheight
    shiftY = handles.imheight-yL;
    borderY = 1:(shiftY+1);
elseif yL == 1
    shiftY = cellsize*2+1-yR;
    borderY = (yL+shiftY):(cellsize*2+1);
end

template(borderY,borderX) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane templateCH -1],tp,handles.channelnames,{[yL yR], [xL xR]},[]);
signalIm(borderY,borderX) = loadimage(handles.filetype,get(handles.edit_fileformat,'String'),[row col field plane nominCH denomCH],tp,handles.channelnames,{[yL yR], [xL xR]},[ThresL ThresH]);

handles.ind_cellpath = ind_cellpath;
handles.template(:,:,tp) = template;
handles.mysignal(:,:,tp) = signalIm;
handles.corners{tp} = [xL xR yL yR];

guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);
updatexyplots(handles);
updatesignalplot(handles);


function edit_shiftsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_shiftsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_shiftsize as text
%        str2double(get(hObject,'String')) returns contents of edit_shiftsize as a double


% --- Executes during object creation, after setting all properties.
function edit_shiftsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_shiftsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fps_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fps as text
%        str2double(get(hObject,'String')) returns contents of edit_fps as a double


% --- Executes during object creation, after setting all properties.
function edit_fps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_pausesignal.
function pushbutton_pausesignal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pausesignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.greenflag = 0;
guidata(hObject, handles);


% --- Executes on button press in pushbutton_playsignal.
function pushbutton_playsignal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_playsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currenttp = str2num(get(handles.edit_signalcurrentframe,'String'));
lasttp = str2num(get(handles.edit_lasttp,'String'));

handles.greenflag = 1;
guidata(hObject, handles);
handles = guidata(hObject);
for tp = currenttp:lasttp
    set(handles.edit_signalcurrentframe,'String',num2str(tp));
    set(handles.slider_signal,'Value',tp);
    updatesignal(handles);
    pause(1/str2num(get(handles.edit_fps,'String')));
    handles = guidata(hObject);
    if handles.greenflag==0
        break;
    end
end
handles.greenflag = 1;
guidata(hObject, handles);

% --- Executes on button press in pushbutton_lastframesignal.
function pushbutton_lastframesignal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_lastframesignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lasttp = str2num(get(handles.edit_lasttp,'String'));

set(handles.edit_signalcurrentframe,'String',num2str(lasttp));
set(handles.slider_signal,'Value',lasttp);
updatesignal(handles);

% --- Executes on button press in pushbutton_firstframesignal.
function pushbutton_firstframesignal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_firstframesignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
firsttp = str2num(get(handles.edit_firsttp,'String'));

set(handles.edit_signalcurrentframe,'String',num2str(firsttp));
set(handles.slider_signal,'Value',firsttp);
updatesignal(handles);


% --- Executes on button press in pushbutton_adoptcell.
function pushbutton_adoptcell_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_adoptcell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sourceF = str2num(get(handles.edit_sourceframe,'String'));
targetF = str2num(get(handles.edit_targetframes,'String'));


cellmask = handles.currentcellmask;

if length(sourceF)==1 & length(targetF)>0

    for i=targetF
        cellmask(:,:,i) = cellmask(:,:,sourceF);
    end

    handles.currentcellmask = cellmask;

    guidata(hObject, handles);
    handles = guidata(hObject);
    updatexyplots(handles);
    updatesignalplot(handles);
    updatetemplate(handles);
    updatesignal(handles);
else
    display('Size of source and/or target frames are not appropriate.');
    return;
end


% --- Executes on button press in togglebutton_selected.
function togglebutton_selected_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_selected
cellNo = str2num(get(handles.edit_cellno,'String'));
myvalue = get(hObject,'Value'); 
selected_cells = handles.selected_cells;
switch myvalue
    case 1
        if isempty(selected_cells)
            selected_cells(1) = cellNo;
        else
            if isempty(find( cellNo == selected_cells ))
                selected_cells(end+1) = cellNo;
                selected_cells = sort(selected_cells);
            end
        end
        
    case 0
       
        foundInd = find( cellNo == selected_cells );
        if ~isempty(foundInd)
            if length(selected_cells)==1
                selected_cells = [];
            elseif foundInd == length(selected_cells)
                selected_cells = selected_cells(1:(foundInd-1));
            elseif foundInd == 1 
                selected_cells = selected_cells(2:end);  
            elseif foundInd>1  & foundInd < length(selected_cells)
                selected_cells = [selected_cells(1:(foundInd-1)) selected_cells((foundInd+1):end)]
            end
        end
end
handles.selected_cells = selected_cells;
guidata(hObject, handles);


% --- Executes on button press in pushbutton_definecyto.
function pushbutton_definecyto_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_definecyto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tp = str2num(get(handles.edit_templatecurrentframe,'String'));
h = imfreehand(handles.axes_template);
accepted_pos = wait(h);
BW = createMask(h);
delete(h);
handles.currentcytomask(:,:,tp) = BW;
guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);
updatesignalplot(handles);


% --- Executes on button press in pushbutton_calcyto.
function pushbutton_calcyto_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calcyto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nucmask  = handles.currentnucmask;
cellmask = handles.currentcellmask;

for tp = str2num(get(handles.edit_cytoframes,'String'));
    handles.currentcytomask(:,:,tp) = bwmorph(cellmask(:,:,tp),'erode',1) & ~ nucmask(:,:,tp);
end
guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);
updatesignalplot(handles);
updatesignalplot(handles);


% --- Executes on button press in pushbutton_adoptcyto.
function pushbutton_adoptcyto_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_adoptcyto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sourceF = str2num(get(handles.edit_sourceframe,'String'));
targetF = str2num(get(handles.edit_targetframes,'String'));


cytomask = handles.currentcytomask;

if length(sourceF)==1 & length(targetF)>0

    for i=targetF
        cytomask(:,:,i) = cytomask(:,:,sourceF);
    end

    handles.currentcytomask = cytomask;

    guidata(hObject, handles);
    handles = guidata(hObject);
    updatexyplots(handles);
    updatesignalplot(handles);
    updatetemplate(handles);
    updatesignal(handles);
else
    display('Size of source and/or target frames are not appropriate.');
    return;
end


% --- Executes on button press in pushbutton_clearnuc.
function pushbutton_clearnuc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clearnuc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for tp = str2num(get(handles.edit_clearframes,'String'));
    BW = zeros(size(handles.currentnucmask,1),size(handles.currentnucmask,2));
    handles.currentnucmask(:,:,tp) = BW;
end
guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);
updatesignalplot(handles);

% --- Executes on button press in pushbutton_clearcell.
function pushbutton_clearcell_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clearcell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for tp = str2num(get(handles.edit_clearframes,'String'));
    BW = zeros(size(handles.currentcellmask,1),size(handles.currentcellmask,2));
    handles.currentcellmask(:,:,tp) = BW;
end
guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);
updatesignalplot(handles);

% --- Executes on button press in pushbutton_clearcyto.
function pushbutton_clearcyto_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clearcyto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for tp = str2num(get(handles.edit_clearframes,'String'));
    BW = zeros(size(handles.currentcytomask,1),size(handles.currentcytomask,2));
    handles.currentcytomask(:,:,tp) = BW;
end
guidata(hObject, handles);
handles = guidata(hObject);
updatetemplate(handles);
updatesignal(handles);
updatesignalplot(handles);


function edit_clearframes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_clearframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_clearframes as text
%        str2double(get(hObject,'String')) returns contents of edit_clearframes as a double


% --- Executes during object creation, after setting all properties.
function edit_clearframes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_clearframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cytoframes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cytoframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cytoframes as text
%        str2double(get(hObject,'String')) returns contents of edit_cytoframes as a double


% --- Executes during object creation, after setting all properties.
function edit_cytoframes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cytoframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_exportdata.
function pushbutton_exportdata_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_exportdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nucmask  = handles.currentnucmask;
cellmask = handles.currentcellmask;
cytomask = handles.currentcytomask;

signalIm = handles.mysignal;
firsttp = str2num(get(handles.edit_firsttp,'String'));
lasttp = str2num(get(handles.edit_lasttp,'String'));

sindex=1;
plotsignal = [];
for tp = firsttp:lasttp
    
    switch handles.segmenttype
        case 1
            BW = nucmask(:,:,tp);
            [Rs Cs] = find(BW);
            if ~isempty(find(BW))
                finalSignal = signalIm(Rs,Cs,tp);
                
                switch handles.stattype
                    case 1
                        plotsignal(sindex,:) = [tp mean(finalSignal(:))];
                    case 2
                        plotsignal(sindex,:) = [tp median(finalSignal(:))];
                    case 3
                        plotsignal(sindex,:) = [tp sum(finalSignal(:))];
                end
                sindex = sindex+1;
            end
            
        case 2
            if ~isempty(find(cytomask(:,:,tp)))
                BW = cytomask(:,:,tp);
            else
                BW = cellmask(:,:,tp) & ~ nucmask(:,:,tp);
            end
            [Rs Cs] = find(BW);
            if ~isempty(find(BW))
                finalSignal = signalIm(Rs,Cs,tp);
                
                switch handles.stattype
                    case 1
                        plotsignal(sindex,:) = [tp mean(finalSignal(:))];
                    case 2
                        plotsignal(sindex,:) = [tp median(finalSignal(:))];
                    case 3
                        plotsignal(sindex,:) = [tp sum(finalSignal(:))];
                end
                sindex = sindex+1;
            end
        case 3
            BW1 = nucmask(:,:,tp);
            
            if ~isempty(find(cytomask(:,:,tp)))
                BW2 = cytomask(:,:,tp);
            else
                BW2 = cellmask(:,:,tp) & ~ nucmask(:,:,tp);
            end
            
            [Rs1 Cs1] = find(BW1);
            [Rs2 Cs2] = find(BW2);
            if ~isempty(find(BW1)) & ~isempty(find(BW2))
                nuclear = signalIm(Rs1,Cs1,tp);
                cytosol = signalIm(Rs2,Cs2,tp);
                
                switch handles.stattype
                    case 1
                        plotsignal(sindex,:) = [tp mean(nuclear(:))/mean(cytosol(:))];
                    case 2
                        plotsignal(sindex,:) = [tp median(nuclear(:))/median(cytosol(:))];
                    case 3
                        plotsignal(sindex,:) = [tp sum(nuclear(:))/sum(cytosol(:))];
                end
                
                sindex = sindex+1;
            end
    end
end

uisave('plotsignal','mydata');



function edit_ndfilename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ndfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ndfilename as text
%        str2double(get(hObject,'String')) returns contents of edit_ndfilename as a double


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
[filename,pathname,FilterIndex] = uigetfile('*.nd', 'Choose metamorph ND file');
if FilterIndex~=0
    set(handles.edit_ndfilename,'String',filename);
    handles.ndfilename = filename;
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

% --- Executes on button press in pushbutton_loadbyndfile.
function pushbutton_loadbyndfile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadbyndfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[notp stagePos stageName waveName] = readndfile(handles.ndfilename);

if notp==-1
    return;
end

set(handles.edit_firsttp,'String',num2str(1));
set(handles.edit_lasttp,'String',num2str(notp));
set(handles.edit_templatecurrentframe,'String',num2str(1));
set(handles.edit_signalcurrentframe,'String',num2str(1));

set(handles.popupmenu_stagePos,'String',stagePos);
set(handles.popupmenu_stagePos,'Value',1);
set(handles.edit_stageInfo,'String',stageName{1});
handles.stageName = stageName;
handles.channelnames = waveName;

prefix = handles.ndfilename(1:(end-3));
handles.prefix = prefix;
fileformat = [prefix '_%s_s1_t%g.tif'];
set(handles.edit_fileformat,'String',fileformat);

tokens   = regexp(stageName{1}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
if length(tokens)>0
    
    row = tokens{1}{1};
    col = tokens{1}{2};
    set(handles.edit_row,'String',row);
    set(handles.edit_col,'String',col);
else
    set(handles.edit_row,'String','1');
    set(handles.edit_col,'String','1');
end

handles.filetype = 3;
set(handles.radiobutton_customim,'Value',1);
set(handles.togglebutton_nucleus,'Value',1);
set(handles.togglebutton_mean,'Value',1);
handles.cellrecogmethod = 1;
handles.segmenttype = 1;
handles.stattype = 1;
handles.greenflag=1;

minF = str2num(get(handles.edit_firsttp,'String'));
maxF = str2num(get(handles.edit_lasttp,'String'));

set(handles.slider_template,'Max',maxF);
set(handles.slider_template,'Min',minF);
set(handles.slider_template,'Value',str2num(get(handles.edit_templatecurrentframe,'String')));
set(handles.slider_template,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);

set(handles.slider_signal,'Max',maxF);
set(handles.slider_signal,'Min',minF);
set(handles.slider_signal,'Value',str2num(get(handles.edit_signalcurrentframe,'String')));
set(handles.slider_signal,'SliderStep',[1/(maxF-minF) 1/(maxF-minF)]);
set(handles.edit_sourceframe,'String',num2str(minF) );
set(handles.edit_targetframes,'String',[num2str(minF) ':' num2str(maxF)]);


guidata(hObject, handles);  

% --- Executes on selection change in popupmenu_stagePos.
function popupmenu_stagePos_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_stagePos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_stagePos contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_stagePos

set(handles.edit_stageInfo,'String',handles.stageName{get(hObject,'Value')});
fileformat = [handles.prefix '_%s_s' num2str(get(hObject,'Value')) '_t%g.tif'];
set(handles.edit_fileformat,'String',fileformat);

tokens   = regexp(handles.stageName{get(hObject,'Value')}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
row = tokens{1}{1};
col = tokens{1}{2};
set(handles.edit_row,'String',row);
set(handles.edit_col,'String',col);


tokens   = regexp(handles.stageName{get(hObject,'Value')}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
if length(tokens)>0
    row = tokens{1}{1};
    col = tokens{1}{2};
    set(handles.edit_row,'String',row);
    set(handles.edit_col,'String',col);
else
    set(handles.edit_row,'String',num2str(get(hObject,'Value')));
    set(handles.edit_col,'String','1');
end

guidata(hObject, handles);  

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
        set(handles.edit_fileformat,'String','r%02.0fc%02.0ff%02.0fp%02.0frc%1.0f-ch1sk%ufk1fl1.tiff');
        handles.filetype = 1;
    case 'radiobutton_tiffstack'
        handles.filetype = 2;
    case 'radiobutton_customim'      
        handles.filetype = 3;
end
guidata(hObject, handles);   


% --- Executes during object creation, after setting all properties.
function togglebutton_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to togglebutton_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit_gateL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gateL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gateL as text
%        str2double(get(hObject,'String')) returns contents of edit_gateL as a double


% --- Executes during object creation, after setting all properties.
function edit_gateL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gateL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gateH_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gateH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gateH as text
%        str2double(get(hObject,'String')) returns contents of edit_gateH as a double


% --- Executes during object creation, after setting all properties.
function edit_gateH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gateH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
