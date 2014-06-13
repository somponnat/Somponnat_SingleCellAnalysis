function varargout = FRETPlotting(varargin)
% Edit the above text to modify the response to help mytab

% Last Modified by GUIDE v2.5 07-Apr-2014 13:45:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FRETPlotting_OpeningFcn, ...
    'gui_OutputFcn',  @FRETPlotting_OutputFcn, ...
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

function outputim = loadsignalV2(handles,channel,tp,signalformat,blankformat)
totalCHs = handles.totalCHs;
channelnames = handles.channelnames;
filetype = handles.filetype;
outputim = [];

if handles.useblank
    switch filetype
        case 1
            
            signal_filename = sprintf(signalformat,channel,tp);
            blank_filename = sprintf(blankformat,channel,tp);
            
            if exist(fullfile(handles.ndpathname,signal_filename),'file')
                signalim = im2double(imread(fullfile(handles.ndpathname,signal_filename)));
                if exist(fullfile(handles.ndpathname,blank_filename),'file')
                    blankim = im2double(imread(fullfile(handles.ndpathname,blank_filename)));
                    smoothed = fftolamopt2(blankim,handles.gaussian,handles.smooth_opt,'same');
                    normim =(smoothed./mean(smoothed(:)));
                    clear blankim smoothed;
                    outputim = signalim./normim;
                else
                    outputim = signalim;
                end
            end
            
            
        case 2
            if exist(fullfile(handles.ndpathname,signalformat),'file')
                outputim = im2double(imread(fullfile(handles.ndpathname,signalformat),'Index',totalCHs*(tp-1)+channel));
            else
                set(handles.edit_commu,'String',[signalformat ' does not exist.']);
            end
        case 3
            signal_filename = sprintf(signalformat,channelnames{channel},tp);
            blank_filename = sprintf(blankformat,channelnames{channel},tp);
            if exist(fullfile(handles.ndpathname,signal_filename),'file')
                signalim = im2double(imread(fullfile(handles.ndpathname,signal_filename)));
                if exist(fullfile(handles.ndpathname,blank_filename),'file')
                    blankim = im2double(imread(fullfile(handles.ndpathname,blank_filename)));
                    smoothed = fftolamopt2(blankim,handles.gaussian,handles.smooth_opt,'same');
                    normim =(smoothed./mean(smoothed(:)));
                    clear blankim smoothed;
                    outputim = signalim./normim;
                else
                    outputim = signalim;
                end
            end
            
    end
else
    
    switch filetype
        case 1
            
            filename = sprintf(signalformat,channel,tp);
            if exist(fullfile(handles.ndpathname,filename),'file')
                outputim = im2double(imread(fullfile(handles.ndpathname,filename)));
            end
        case 2
            if exist(signalformat,'file')
                outputim = im2double(imread(fullfile(handles.ndpathname,signalformat),'Index',totalCHs*(tp-1)+channel));
            end
        case 3
            filename = sprintf(signalformat,channelnames{channel},tp);
            if exist(fullfile(handles.ndpathname,filename),'file');
                outputim = im2double(imread(fullfile(handles.ndpathname,filename)));
            end
    end
end

function outputim = loadblankV2(handles,channel,tp,blankformat)
filetype = handles.filetype;
channelnames = handles.channelnames;
totalCHs = handles.totalCHs;
outputim = [];

switch filetype
    case 1
        
        filename = sprintf(blankformat,channel,tp);
        if exist(fullfile(handles.ndpathname,filename),'file')
            outputim = im2double(imread(fullfile(handles.ndpathname,filename)));
        else
            set(handles.edit_commu,'String',[filename ' does not exist.']);
        end
    case 2
        if exist(signalformat,'file')
            outputim = im2double(imread(fullfile(handles.ndpathname,blankformat),'Index',totalCHs*(tp-1)+channel));
        else
            set(handles.edit_commu,'String',[signalformat ' does not exist.']);
        end
    case 3
        filename = sprintf(blankformat,channelnames{channel},tp);
        if exist(fullfile(handles.ndpathname,filename),'file');
            outputim = im2double(imread(fullfile(handles.ndpathname,filename)));
        else
            set(handles.edit_commu,'String',[filename ' does not exist.']);
        end
end

function ratioIm = calculateFRETV2(handles,bg,tp,nominCH,denominCH,signalformat,blankformat)

filterParam1 = handles.filterParam1;
filterParam2 = handles.filterParam2;
bgsize = handles.bgsize;
signalShiftN = handles.signalShiftN;
signalShiftD = handles.signalShiftD;

nominIM = loadsignalV2(handles,nominCH,tp,signalformat,blankformat);
if denominCH == -1
    denomIM = im2double(ones(size(nominIM)));
else
    denomIM = loadsignalV2(handles,denominCH,tp,signalformat,blankformat);
end

% Load blank images if user chose to use BG points from BLANK images
if handles.bgnominType==3
    nominBLK = loadblankV2(handles,nominCH,tp,blankformat);
    
end
if handles.bgdenominType==3
    if denominCH == -1
        denomBLK = im2double(ones(size(nominIM)));
    else
        denomBLK = loadblankV2(handles,denominCH,tp,blankformat);
    end
end


if handles.illumlogic
    normN = ifft2(ifftshift(fftshift(fft2(nominIM)).*hbutter(nominIM,filterParam1,filterParam2)));
    normD = ifft2(ifftshift(fftshift(fft2(denomIM)).*hbutter(denomIM,filterParam1,filterParam2)));
else
    normN = nominIM;
    normD = denomIM;
end
BG_N = 0;
BG_D = 0;

switch handles.bgnominType
    case 2
        if ~isempty(bg)
            for b=1:size(bg{tp},1)
                xL=max(bg{tp}(b,1)-bgsize,1);
                xR=min(bg{tp}(b,1)+bgsize,size(nominIM,2));
                yL=max(bg{tp}(b,2)-bgsize,1);
                yR=min(bg{tp}(b,2)+bgsize,size(nominIM,1));
                selectedN = normN(yL:yR,xL:xR);
                BG_N(b) = mean(selectedN(:));
            end
        end
    case 3
        xL=round(size(nominIM,2)/2)-2*bgsize;
        xR=round(size(nominIM,2)/2)+2*bgsize;
        yL=round(size(nominIM,1)/2)-2*bgsize;
        yR=round(size(nominIM,1)/2)+2*bgsize;
        BG_N = mean(nominBLK(:));
end

switch handles.bgdenominType
    case 2
        if ~isempty(bg)
            for b=1:size(bg{tp},1)
                xL=max(bg{tp}(b,1)-bgsize,1);
                xR=min(bg{tp}(b,1)+bgsize,size(nominIM,2));
                yL=max(bg{tp}(b,2)-bgsize,1);
                yR=min(bg{tp}(b,2)+bgsize,size(nominIM,1));
                selectedD = normD(yL:yR,xL:xR);
                BG_D(b) = mean(selectedD(:));
            end
        end
    case 3
        xL=round(size(nominIM,2)/2)-2*bgsize;
        xR=round(size(nominIM,2)/2)+2*bgsize;
        yL=round(size(nominIM,1)/2)-2*bgsize;
        yR=round(size(nominIM,1)/2)+2*bgsize;
        BG_D = mean(denomBLK(:));
end

normN = normN-mean(BG_N);
normD = normD-mean(BG_D);

normN = normN+signalShiftN;
normD = normD+signalShiftD;

ratioIm =  normN./normD;

function outputim = loadsignal(handles,channel,tp)

filetype = handles.filetype;
signalformat = get(handles.edit_signalformat,'String');
blankformat = get(handles.edit_fileformatBG,'String');
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channelnames = handles.channelnames;
totalCHs = str2num(get(handles.edit_totalCHs,'String'));
outputim = [];

if handles.useblank
    switch filetype
        case 1
            
            signal_filename = sprintf(signalformat,channel,tp);
            blank_filename = sprintf(blankformat,channel,tp);
            
            if exist(fullfile(handles.ndpathname,signal_filename),'file')
                signalim = im2double(imread(fullfile(handles.ndpathname,signal_filename)));
                if exist(fullfile(handles.ndpathname,blank_filename),'file')
                    blankim = im2double(imread(fullfile(handles.ndpathname,blank_filename)));
                    smoothed = fftolamopt2(blankim,handles.gaussian,handles.smooth_opt,'same');
                    normim =(smoothed./mean(smoothed(:)));
                    clear blankim smoothed;
                    outputim = signalim./normim;
                else
                    set(handles.edit_commu,'String',[blank_filename ' does not exist. Use SignalIM only.']);
                    outputim = signalim;
                end
            else
                set(handles.edit_commu,'String',[signal_filename ' does not exist.']);
            end
            
            
        case 2
            set(handles.edit_commu,'String','This function is not available for tiff stack input.');
            if exist(fullfile(handles.ndpathname,signalformat),'file')
                outputim = im2double(imread(fullfile(handles.ndpathname,signalformat),'Index',totalCHs*(tp-1)+channel));
            else
                set(handles.edit_commu,'String',[signalformat ' does not exist.']);
            end
        case 3
            signal_filename = sprintf(signalformat,channelnames{channel},tp);
            blank_filename = sprintf(blankformat,channelnames{channel},tp);
            if exist(fullfile(handles.ndpathname,signal_filename),'file')
                signalim = im2double(imread(fullfile(handles.ndpathname,signal_filename)));
                if exist(fullfile(handles.ndpathname,blank_filename),'file')
                    blankim = im2double(imread(fullfile(handles.ndpathname,blank_filename)));
                    smoothed = fftolamopt2(blankim,handles.gaussian,handles.smooth_opt,'same');
                    normim =(smoothed./mean(smoothed(:)));
                    clear blankim smoothed;
                    outputim = signalim./normim;
                else
                    set(handles.edit_commu,'String',[blank_filename ' does not exist. Use SignalIM only.']);
                    outputim = signalim;
                end
            else
                set(handles.edit_commu,'String',[signal_filename ' does not exist.']);
            end
            
    end
else
    
    switch filetype
        case 1
            
            filename = sprintf(signalformat,channel,tp);
            if exist(fullfile(handles.ndpathname,filename),'file')
                outputim = im2double(imread(fullfile(handles.ndpathname,filename)));
            else
                set(handles.edit_commu,'String',[filename ' does not exist.']);
            end
        case 2
            if exist(signalformat,'file')
                outputim = im2double(imread(fullfile(handles.ndpathname,signalformat),'Index',totalCHs*(tp-1)+channel));
            else
                set(handles.edit_commu,'String',[signalformat ' does not exist.']);
            end
        case 3
            filename = sprintf(signalformat,channelnames{channel},tp);
            if exist(fullfile(handles.ndpathname,filename),'file');
                outputim = im2double(imread(fullfile(handles.ndpathname,filename)));
            else
                set(handles.edit_commu,'String',[filename ' does not exist.']);
            end
    end
end



function outputim = loadblank(handles,channel,tp)
filetype = handles.filetype;
blankformat = get(handles.edit_fileformatBG,'String');
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channelnames = handles.channelnames;
totalCHs = str2num(get(handles.edit_totalCHs,'String'));
outputim = [];

switch filetype
    case 1
        
        filename = sprintf(blankformat,channel,tp);
        if exist(fullfile(handles.ndpathname,filename),'file')
            outputim = im2double(imread(fullfile(handles.ndpathname,filename)));
        else
            set(handles.edit_commu,'String',[filename ' does not exist.']);
        end
    case 2
        if exist(signalformat,'file')
            outputim = im2double(imread(fullfile(handles.ndpathname,blankformat),'Index',totalCHs*(tp-1)+channel));
        else
            set(handles.edit_commu,'String',[signalformat ' does not exist.']);
        end
    case 3
        filename = sprintf(blankformat,channelnames{channel},tp);
        if exist(fullfile(handles.ndpathname,filename),'file');
            outputim = im2double(imread(fullfile(handles.ndpathname,filename)));
        else
            set(handles.edit_commu,'String',[filename ' does not exist.']);
        end
end


% --- Executes just before mytab is made visible.
function FRETPlotting_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mytab (see VARARGIN)

% Choose default command line output for mytab
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
            trackinginfo = varargin{2};  % [nucCH tp_1 tp_end]
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
    set(handles.edit_commu,'String','Need to provide type (input#1) of input file(s)');
    return;
end



switch filetype
    case 1 % PE tiff input
        
        if isempty(imagelocation)
            set(handles.edit_commu,'String','Required image location (input#3) for PE-tiff format');
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
        
        nucCH = trackinginfo(1);
        set(handles.edit_nucCH,'String',num2str(nucCH));
        set(handles.edit_firstframe,'String',num2str(trackinginfo(2)));
        set(handles.edit_lastframe,'String',num2str(trackinginfo(3)));
        set(handles.edit_currentFrame,'String',num2str(trackinginfo(2)));
        
        
        handles.channelnames = channelnames;
        set(handles.edit_row,'String',num2str(row));
        set(handles.edit_col,'String',num2str(col));
        set(handles.edit_field,'String',num2str(field));
        set(handles.edit_plane,'String',num2str(plane));
        set(handles.edit_signalformat,'String',fileformat);
        
    case 2 % Tiff stack input
        
        if isempty(channelnames)
            display('Required channelnames (input#4) for tiffstack input');
            return;
        end
        
        if isempty(fileformat)
            set(handles.edit_commu,'String','Required filename (input#5) for tiffstack input');
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
        
        
        nucCH = trackinginfo(1);
        set(handles.edit_nucCH,'String',num2str(nucCH));
        set(handles.edit_firstframe,'String',num2str(trackinginfo(2)));
        set(handles.edit_lastframe,'String',num2str(trackinginfo(3)));
        set(handles.edit_currentFrame,'String',num2str(trackinginfo(2)));
        
        handles.channelnames = channelnames;
        set(handles.edit_row,'String',num2str(row));
        set(handles.edit_col,'String',num2str(col));
        set(handles.edit_field,'String',num2str(field));
        set(handles.edit_plane,'String',num2str(plane));
        set(handles.edit_signalformat,'String',fileformat);
    case 3 % tif file with custom naming input
        if isempty(fileformat)
            set(handles.edit_commu,'String','Required fileformat (input#5) for custom image input');
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
            [notp stagePos stageName channelnames] = readndfile(pwd,fileformat);
            handles.ndfilename = fileformat;
            handles.ndpathname = pwd;
            if notp==-1
                handles.initialframe = [];
                set(handles.edit_commu,'String',[fileformat ' does not exist. Please re-define input ND file.']);
                guidata(hObject, handles);
                return;
            end
            
            set(handles.edit_firstframe,'String',num2str(1));
            set(handles.edit_lastframe,'String',num2str(notp));
            set(handles.edit_currentFrame,'String',num2str(1));
            nucCH = 1;
            set(handles.popupmenu_stagePosData,'String',stagePos);
            set(handles.popupmenu_stagePosData,'Value',1);
            set(handles.popupmenu_stagePosBG,'String',stagePos);
            set(handles.popupmenu_stagePosBG,'Value',1);
            set(handles.edit_stageInfodata,'String',stageName{1});
            set(handles.edit_stageInfobg,'String',stageName{1});
            handles.stageName = stageName;
            handles.channelnames = channelnames;
            
            prefix = fileformat(1:(end-3));
            handles.prefix = prefix;
            fileformat = [prefix '_%s_s1_t%g.tif'];
            set(handles.edit_signalformat,'String',fileformat);
            
            fileformat = [prefix '_%s_s1_t%g.TIF'];
            set(handles.edit_fileformatBG,'String',fileformat);
            
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
            set(handles.edit_signalformat,'String',fileformat);
            nucCH = trackinginfo(1);
            set(handles.edit_nucCH,'String',num2str(nucCH));
            set(handles.edit_firstframe,'String',num2str(trackinginfo(2)));
            set(handles.edit_lastframe,'String',num2str(trackinginfo(3)));
            set(handles.edit_currentFrame,'String',num2str(trackinginfo(2)));
        end
        
end

handles.filetype = filetype;
handles.framestamp = 1;
set(handles.radiobutton_frameno,'Value',1);
set(handles.listbox_cells,'String',[]);
handles.cellpath = [];
handles.bg = [];
handles.sisterList = [];
handles.greenflag=1;
handles.nucmask = [];
handles.cytomask = [];
handles.cellmask = [];
handles.selected_cells = [];
handles.maskpathname = [];
handles.celltrackpathname = [];
set(handles.edit_mathL,'String',num2str(0.7));
set(handles.edit_mathH,'String',num2str(1.7));
set(handles.edit_templateL,'String',num2str(0));
set(handles.edit_templateH,'String',num2str(1));
set(handles.edit_ch1L,'String',num2str(0));
set(handles.edit_ch1H,'String',num2str(1));
set(handles.edit_ch2L,'String',num2str(0));
set(handles.edit_ch2H,'String',num2str(1));
set(handles.edit_ch3L,'String',num2str(0));
set(handles.edit_ch3H,'String',num2str(1));
handles.celltrackpathname = [];
handles.fcn1 = makeConstrainToRectFcn('impoint',get(handles.axes1,'XLim'),get(handles.axes1,'YLim'));
handles.useblank = 0;
handles.gaussian = 1/273*[1 4 7 4 1;4 16 26 26 4;7 26 41 26 7;4 16 26 16 4;1 4 7 4 1];

guidata(hObject, handles);
handles = guidata(hObject);

fftw('planner', 'hybrid');
load MyColormaps;
handles.mycmap1 = mycmap1;
handles.mycmap2 = mycmap2;
handles.mycmap3 = mycmap3;
handles.mycmap4 = mycmap4;

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));

CH1= handles.CH1;
CH2= handles.CH2;
CH3 = handles.CH3;
cellsize = str2num(get(handles.edit_cellsize,'String'));
tp=str2num(get(handles.edit_currentFrame,'String'));

switch get(handles.popupmenu_nomin,'Value')
    case 1
        nominCH=CH1;
    case 2
        nominCH=CH2;
    case 3
        nominCH=CH3;
end

switch get(handles.popupmenu_denomin,'Value')
    case 1
        denominCH=CH1;
    case 2
        denominCH=CH2;
    case 3
        denominCH=CH3;
end
startImageIndex = 1;
handles.ImageIndex = startImageIndex;
handles.overlayIndex = startImageIndex;
if filetype~=0
    
    template = loadsignal(handles,nucCH,tp);
    load fftexecutiontimes;
    handles.smooth_opt = detbestlength2(FFTrv,FFTiv,IFFTiv,size(template),size(handles.gaussian),isreal(template),isreal(handles.gaussian));
    switch startImageIndex
        
        case 3
            displayIM = loadsignal(handles,CH1,tp);
            imshow(displayIM,[],'Parent',handles.axes1);colormap gray;drawnow;
        case 4
            displayIM = loadsignal(handles,CH2,tp);
            imshow(displayIM,[],'Parent',handles.axes1);colormap gray;drawnow;
        case 5
            displayIM = loadsignal(handles,CH3,tp);
            imshow(displayIM,[],'Parent',handles.axes1);colormap gray;drawnow;
        case 2
            displayIM = template;
            imshow(displayIM,[],'Parent',handles.axes1);colormap gray;drawnow;
        case 1
            displayIM = calculateFRET(handles,tp,nominCH,denominCH,[]);
            displayratioIm(displayIM,handles);
    end
    drawnow;
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
nucCH= str2num(get(handles.edit_nucCH,'String'));
CH1= handles.CH1;
CH2= handles.CH2;
CH3 = handles.CH3;
cellsize = str2num(get(handles.edit_cellsize,'String')); %#ok<*NASGU>
c_tp=str2num(get(handles.edit_firstframe,'String')); %#ok<*ST2NM>
set(handles.edit_currentFrame,'String',num2str(c_tp));
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

totalCHs = str2num(get(handles.edit_totalCHs,'String'));
fileformat = get(handles.edit_signalformat,'String');



switch get(handles.popupmenu_nomin,'Value')
    case 1
        nominCH=CH1;
    case 2
        nominCH=CH2;
    case 3
        nominCH=CH3;
end

switch get(handles.popupmenu_denomin,'Value')
    case 1
        denominCH=CH1;
    case 2
        denominCH=CH2;
    case 3
        denominCH=CH3;
    case 4
        denominCH=-1;
end
handles.gaussian = 1/273*[1 4 7 4 1;4 16 26 26 4;7 26 41 26 7;4 16 26 16 4;1 4 7 4 1];
guidata(hObject, handles);
handles = guidata(hObject);
template = loadsignal(handles,nucCH,c_tp);
load fftexecutiontimes;
handles.smooth_opt = detbestlength2(FFTrv,FFTiv,IFFTiv,size(template),size(handles.gaussian),isreal(template),isreal(handles.gaussian));
clear FFTrv FFTiv IFFTiv;
guidata(hObject, handles);
handles = guidata(hObject);

startImageIndex = 1;
handles.ImageIndex = startImageIndex;
handles.overlayIndex = startImageIndex;
set(handles.togglebutton_math,'Value',1);
set(handles.radiobutton_mathoverlay,'Value',1);
load MyColormaps;
handles.mycmap1 = mycmap1;
handles.mycmap2 = mycmap2;
handles.mycmap3 = mycmap3;
handles.mycmap4 = mycmap4;
handles.outputsignalname = get(handles.edit_outputname,'String');
clear mycmap1 mycmap2 mycmap3 mycmap4
switch startImageIndex
    case 3
        displayIM = loadsignal(handles,CH1,c_tp);
        if ~isempty(displayIM)
            imshow(displayIM,[],'Parent',handles.axes1);colormap gray;drawnow;
        end
    case 4
        displayIM = loadsignal(handles,CH2,c_tp);
        if ~isempty(displayIM)
            imshow(displayIM,[],'Parent',handles.axes1);colormap gray;drawnow;
        end
    case 5
        displayIM = loadsignal(handles,CH3,c_tp);
        if ~isempty(displayIM)
            imshow(displayIM,[],'Parent',handles.axes1);colormap gray;drawnow;
        end
    case 2
        displayIM = template;
        if ~isempty(displayIM)
            imshow(displayIM,[],'Parent',handles.axes1);colormap gray;drawnow;
        end
    case 1
        displayIM = calculateFRET(handles,c_tp,nominCH,denominCH,[]);
        if ~isempty(displayIM)
            displayratioIm(displayIM,handles);
        end
end
set(handles.edit_cellNo,'String',1);
set(handles.listbox_cells,'String',[]);
set(handles.listbox_cells,'Value',1);
handles.cellpath = [];
handles.bg = [];
handles.sisterList = [];

H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
fileattrib(fullfile(handles.ndpathname,H5filename),'+w');
cellpath_name = ['/field' num2str(field) '/cellpath'];
sisterList_name = ['/field' num2str(field) '/sisterList'];
bg_name = ['/field' num2str(field) '/bg'];

if exist(fullfile(handles.ndpathname,H5filename),'file')
    fid = H5F.open(fullfile(handles.ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
        H5F.close(fid);
        cellpathinfo = h5info(fullfile(handles.ndpathname,H5filename), cellpath_name);
        cellpath_mat = h5read(fullfile(handles.ndpathname,H5filename),cellpath_name,[1 1 1], [cellpathinfo.Dataspace.Size(1) cellpathinfo.Dataspace.Size(2) cellpathinfo.Dataspace.Size(3)]);
        
        for tp=first_tp:last_tp
            cellpath{tp} = cellpath_mat(:,:,tp);
        end
        handles.cellpath=cellpath;
        set(handles.edit_commu,'String',[H5filename ' ' cellpath_name 'loaded.']);
        
    else
        handles.cellpath=[];
        set(handles.edit_commu,'String',[H5filename ' ' cellpath_name ' does not exist.']);
    end
    
    
    fid = H5F.open(fullfile(handles.ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if H5L.exists(fid,sisterList_name,'H5P_DEFAULT')
        H5F.close(fid);
        sisterListinfo = h5info(fullfile(handles.ndpathname,H5filename), sisterList_name);
        sisterList_mat = h5read(fullfile(handles.ndpathname,H5filename),sisterList_name,[1 1 1], [sisterListinfo.Dataspace.Size(1) sisterListinfo.Dataspace.Size(2) sisterListinfo.Dataspace.Size(3)]);
        
        if size(sisterList_mat,1) == size(cellpath_mat,1) && size(sisterList_mat,2) == size(cellpath_mat,2) && size(sisterList_mat,3) == size(cellpath_mat,3)
            
            for tp=first_tp:sisterListinfo.Dataspace.Size(3)
                sisterList{tp} = sisterList_mat(:,:,tp);
            end
        else
            for tp=first_tp:cellpathinfo.Dataspace.Size(3)
                sisterList{tp} = -1*ones(size(cellpath_mat,1),3);
            end
        end
        
        handles.sisterList=sisterList;
        
    else
        if length(cellpath{tp}) > 0
            for tp=first_tp:last_tp
                sisterList{tp} = -1*ones(size(cellpath{tp},1),3);
            end
            handles.sisterList=sisterList;
        else
            handles.sisterList = [];
        end
    end
    
    
    fid = H5F.open(fullfile(handles.ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if H5L.exists(fid,bg_name,'H5P_DEFAULT')
        H5F.close(fid);
        bginfo = h5info(fullfile(handles.ndpathname,H5filename), bg_name);
        bg_mat = h5read(fullfile(handles.ndpathname,H5filename),bg_name,[1 1 1], [bginfo.Dataspace.Size(1) bginfo.Dataspace.Size(2) bginfo.Dataspace.Size(3)]);
        if size(bg_mat,3) == size(cellpath_mat,3)
            
            for tp=first_tp:last_tp
                bg{tp} = bg_mat(:,:,tp);
            end
            handles.bg = bg;
        else
            handles.bg=[];
        end
    else
        handles.bg=[];
    end
else
    handles.cellpath=[];
    handles.sisterList = [];
    handles.bg=[];
    set(handles.edit_commu,'String',[H5filename ' does not exist.']);
end

guidata(hObject, handles);
handles = guidata(hObject);
updatecurrentImage(handles,c_tp);
updateOutputList(handles);
handles.greenflag=1;
handles.nucmask = [];
handles.cytomask = [];
handles.cellmask = [];
handles.selected_cells = [];

handles.fcn1 = makeConstrainToRectFcn('impoint',get(handles.axes1,'XLim'),get(handles.axes1,'YLim'));
guidata(hObject, handles);

% --- Executes on button press in pushbutton_editpoint.
function pushbutton_editpoint_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_editpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellpath = handles.cellpath;
sisterList = handles.sisterList;
coord = get(handles.edit_coord,'String');

tp = str2num(get(handles.edit_currentFrame,'String'));
cell   = str2num(get(handles.edit_cellNo,'String'));
sister = str2num(get(handles.edit_sister,'String'));
x = str2num(regexp(coord,'\d+(?=,\d+)','match'));
y = str2num(regexp(coord,'(?<=\d+,)\d+','match'));

cellpath{tp}(cell,:) = [x y];
sisterList{tp}(cell,:) = sister;
handles.sisterList = sisterList;
handles.cellpath = cellpath;
guidata(hObject, handles);

function [x y nucMask cytoMask] = templateToCentroid(handles,M,xg,yg,mycytosize)

maxI = handles.maxI;
invertLog = handles.invertedLog;
outerbox = str2num(get(handles.edit_cellsize,'String'));

if invertLog
    inverted = ((maxI)-1)-M;
end
BWc = zeros(size(M));
for i=2:0.5:6
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


nucMask = bwselect(BW,xg,yg);
BW3 = bwmorph(nucMask,'dilate',mycytosize);
cytoMask = BW3-nucMask;

function [x y nucMask cytoMask] = templateToCentroid2(handles,M,xg,yg,mycytosize)

maxI = handles.maxI;
invertLog = handles.invertedLog;
outerbox = str2num(get(handles.edit_cellsize,'String'));

if invertLog
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
S  = regionprops(outBW, 'centroid');

if ~isfield(S,'Centroid') || length(S)>1
    x = xg;
    y = yg;
    
else
    x = round(S.Centroid(1));
    y = round(S.Centroid(2));
end

nucMask = bwselect(BW,xg,yg);
BW3 = bwmorph(nucMask,'dilate',mycytosize);
cytoMask = BW3-nucMask;


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
    old_{i} = old(:,:,i); %#ok<*AGROW>
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
    if Q1<=dt*.18^2 && Q2<=dt*.18^2 %#ok<*BDSCI>
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


% --- Executes on button press in pushbutton_corrsetting.
function pushbutton_corrsetting_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_corrsetting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    cellpath = handles.cellpath;
    cellsize = str2num(get(handles.edit_cellsize,'String'));
    tp = str2num(get(handles.edit_currentFrame,'String'));
    cell = get(handles.listbox_cells,'Value');
    fcn = makeConstrainToRectFcn('imrect',get(handles.axes1,'XLim'),get(handles.axes1,'YLim'));
    h1 = imrect(handles.axes1, [cellpath{tp}(cell,1)-cellsize cellpath{tp}(cell,2)-cellsize 2*cellsize-1 2*cellsize-1]);
    setPositionConstraintFcn(h1,fcn);
    setResizable(h1,0);
    handles.rect1 = h1;
    
elseif button_state == get(hObject,'Min')
    delete(handles.rect1);
end

guidata(hObject, handles);


% --- Executes on selection change in listbox_cells.
function listbox_cells_Callback(hObject, ~, handles)
% hObject    handle to listbox_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_cells contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_cells
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

selected_cell = get(hObject,'Value');
set(handles.edit_cellNo,'String',num2str(selected_cell));
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
nucCH= str2num(get(handles.edit_nucCH,'String'));
cellpath = handles.cellpath;
sisterList = handles.sisterList;
cellsize = str2num(get(handles.edit_cellsize,'String'));
tp = str2num(get(handles.edit_currentFrame,'String'));

H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];

maskdatasetname = ['/field' num2str(field) '/segmentsCH' num2str(nucCH)];
fid = H5F.open(fullfile(handles.ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,maskdatasetname,'H5P_DEFAULT')
    H5F.close(fid);
    maskinfo = h5info(fullfile(handles.ndpathname,H5filename), maskdatasetname);
    
    selectedcells_name = ['/field' num2str(field) '/selectedcells'];
    selected_cells = h5read(fullfile(handles.ndpathname,H5filename),selectedcells_name);
    
    set(handles.edit_collectedcells,'String',num2str(selected_cells'));

    
    firsttp = str2num(get(handles.edit_firstframe,'String'));
    lasttp = str2num(get(handles.edit_lastframe,'String'));
    currentframe = loadsignal(handles,nucCH,tp);
    imwidth = size(currentframe,2);
    imheight = size(currentframe,1);
    [new_cellpath,new_sisterList] = removeSister(cellpath,sisterList,firsttp,lasttp,1:length(cellpath{lasttp}));
    ind_cellpath = pos_path(new_cellpath,new_sisterList,selected_cell,firsttp,lasttp,imheight,imwidth);
    handles.ind_cellpath = ind_cellpath;
    allmasks =permute(h5read(fullfile(handles.ndpathname,H5filename),maskdatasetname,[selected_cell firsttp 1 1 1], [1 lasttp-firsttp+1 3 maskinfo.Dataspace.Size(4) maskinfo.Dataspace.Size(5)]),[4 5 2 3 1]);
    
    nucmask  = double(allmasks(:,:,:,1));
    cellmask  = double(allmasks(:,:,:,2));
    cytomask  = double(allmasks(:,:,:,3));
    
    if cellsize~=(size(nucmask(:,:,tp),1)-1)/2
        cellsize = (size(nucmask(:,:,tp),1)-1)/2;
    end
    
    if tp > size(ind_cellpath,1)
        set(handles.edit_commu,'String',['Cell#' num2str(selected_cell) ' is dead at frame' num2str(tp)]);
    else
        
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
        
        
        cell_template = zeros(2*cellsize+1,2*cellsize+1);
        cell_template(borderY,borderX) = imadjust(currentframe(yL:yR,xL:xR));
        cp = zeros(size(cell_template));
        cp(cellsize,cellsize) = 1;
        imAdj = cell_template-cp;
        imOut=cat(3,max(imAdj,cp),imAdj,imAdj);
        
        
        if ~isempty(cellpath) %#ok<*USENS>
            nosisterInd = find(sisterList{tp}(:,1)==-1);
            sister1Ind = find(sisterList{tp}(:,1)~=-1);
            updatecurrentImage(handles,tp);
            axes(handles.axes1);hold on;
            plot(ind_cellpath(tp,1),ind_cellpath(tp,2),'ro');
            hold off;
        end
        
        if ~isempty(find(selected_cell==selected_cells', 1)) && ~isempty(find(nucmask(:,:,tp)==1,1))
            set(handles.edit_commu,'String',['Cell#' num2str(selected_cell) ' was selected']);
            NucMask  = nucmask(:,:,tp) ;
            
            if ~isempty(find(cellmask(:,:,tp),1))
                CellMask = cellmask(:,:,tp);
            else
                CellMask = bwmorph(NucMask,'dilate',str2num(get(handles.edit_cytosize,'String')));
                cellmask(:,:,tp) = CellMask;
                
            end
            
            if ~isempty(find(cytomask(:,:,tp),1))
                CytoMask = cytomask(:,:,tp);
            else
                CytoMask = bwmorph(CellMask,'erode',1)-NucMask;
                cytomask(:,:,tp) = CytoMask;
                
            end
            
            bwNucEdge = bwperim(NucMask);
            bwCytoEdge = bwperim(CytoMask);
            bwCellEdge = bwperim(CellMask);
            
            bwFinal = bwNucEdge| bwCellEdge | bwCytoEdge;
            imAdj = imAdj-bwFinal;
            imOut=cat(3,max(imAdj,(bwNucEdge | cp)),max(imAdj,bwCellEdge),max(imAdj,bwCytoEdge));
        else
            set(handles.edit_commu,'String',['Cell#' num2str(selected_cell) ' has no masks']);
        end
        imshow(imOut,[],'Parent',handles.axes2);
        
        set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',ind_cellpath(tp,1),ind_cellpath(tp,2))));
        
        set(handles.edit_sister,'String',num2str(sisterList{tp}(selected_cell,:)));
        
    end
    
    handles.nucmask = nucmask;
    handles.cellmask = cellmask;
    handles.cytomask = cytomask;
    
    guidata(hObject, handles);
else
    H5F.close(fid);
    if ~isempty(cellpath) %#ok<*USENS>
        updatecurrentImage(handles,tp);
        axes(handles.axes1);hold on;
        plot(cellpath{tp}(selected_cell,1),cellpath{tp}(selected_cell,2),'ro');hold off;
        set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',cellpath{tp}(selected_cell,1),cellpath{tp}(selected_cell,2))));
        set(handles.edit_sister,'String',num2str(sisterList{tp}(selected_cell,:)));
        
        set(handles.edit_commu,'String',['Cell#' num2str(selected_cell) ' was selected']);
        nucFolder = [handles.ndpathname filesep 'nuclearMask'];
        cellFolder = [handles.ndpathname filesep 'cellMask'];
        nucCH = str2num(get(handles.edit_nucCH,'String'));
        cellCH = str2num(get(handles.edit_cellCH,'String'));
        signalformat = get(handles.edit_signalformat,'String');
        channelnames = handles.channelnames;
        if exist(fullfile(nucFolder, sprintf(signalformat,channelnames{nucCH},tp)),'file')
            nuc_im   = imread(fullfile(nucFolder, sprintf(signalformat,channelnames{nucCH},tp)));
            if exist(fullfile(cellFolder,sprintf(signalformat,channelnames{cellCH},tp)),'file')
                cell_im  = imread(fullfile(cellFolder,sprintf(signalformat,channelnames{cellCH},tp)));
            elseif exist(fullfile(cellFolder,sprintf(signalformat,channelnames{nucCH},tp)),'file')
                cell_im  = imread(fullfile(cellFolder,sprintf(signalformat,channelnames{nucCH},tp)));
            end
            imwidth = size(nuc_im,2);
            imheight = size(nuc_im,1);
            cX = cellpath{tp}(selected_cell,1);
            cY = cellpath{tp}(selected_cell,2);
            cellmaskT = bwselect(cell_im,cX,cY,8);
            nucmaskT  = bwselect(nuc_im,cX,cY,8);
            if ~isempty(find(nucmaskT,1)) || ~isempty(find(cellmaskT,1))
                S = regionprops(cellmaskT,{'BoundingBox'});
                BoundingBox = round(S(1).BoundingBox);
                
                xL=max(BoundingBox(1),1);
                xR=min(xL+BoundingBox(3),imwidth);
                yL=max(BoundingBox(2),1);
                yR=min(yL+BoundingBox(4),imheight);
                
                cellmask = cellmaskT(yL:yR,xL:xR);
                nucmask  = nucmaskT(yL:yR,xL:xR);
                cytomask = cellmask-nucmask;
                bwNucEdge = double(bwperim(nucmask));
                bwCytoEdge = double(bwperim(cytomask));
                bwCellEdge = double(bwperim(cellmask));
                handles.nucmask = nucmask;
                handles.cellmask = cellmask;
                handles.cytomask = cytomask;
                guidata(hObject, handles);
                imOut=cat(3,bwNucEdge,bwCellEdge,bwCytoEdge);
                imshow(imOut,'Parent',handles.axes2);
            end
            
        end

    end
end


function  out_cellpath = pos_path(cellpath,sisterList,cellNo,firsttp,input_lasttp,imheight,imwidth)
% check if this path is good
%#1 if acquiring and being inferior sister, combine trace of self with
%prior-sister
% has sister?

ind_cellpath = zeros(input_lasttp,2);

for t=firsttp:input_lasttp
    ind_cellpath(t,:) = cellpath{t}(cellNo,:);
end
deathInd = find(ind_cellpath(:,1)==-2,1,'first');

if isempty(deathInd)
    out_cellpath = ind_cellpath;
else
    out_cellpath = ind_cellpath(1:(deathInd-1),:);
end

% --- Outputs from this function are returned to the command line.
function varargout = FRETPlotting_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_tofirstframe.
function pushbutton_tofirstframe_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_tofirstframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

tp=str2num(get(handles.edit_firstframe,'String'));
updatecurrentImage(handles,tp);
guidata(hObject, handles);


% --- Executes on button press in pushbutton_topreviousframe.
function pushbutton_topreviousframe_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_topreviousframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end



c_tp = str2num(get(handles.edit_currentFrame,'String'));
first_tp = str2num(get(handles.edit_firstframe,'String'));
if c_tp-1<first_tp
    tp=first_tp;
else
    tp=c_tp-1;
end

updatecurrentImage(handles,tp);
guidata(hObject, handles);



% --- Executes on button press in pushbutton_tonextframe.
function pushbutton_tonextframe_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_tonextframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

c_tp = str2num(get(handles.edit_currentFrame,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

if c_tp+1>last_tp
    tp=last_tp;
else
    tp=c_tp+1;
end

updatecurrentImage(handles,tp);
guidata(hObject, handles);

% --- Executes on button press in pushbutton_tolastframe.
function pushbutton_tolastframe_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_tolastframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end
tp=str2num(get(handles.edit_lastframe,'String'));
updatecurrentImage(handles,tp);
guidata(hObject, handles);

function edit_currentFrame_Callback(hObject, ~, handles)
% hObject    handle to edit_currentFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_currentFrame as text
%        str2double(get(hObject,'String')) returns contents of edit_currentFrame as a double
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end
tp=round(str2num(get(hObject,'String')));
first_tp = str2num(get(handles.edit_firstframe,'String'));
updatecurrentImage(handles,tp);
guidata(hObject, handles);

function [forwardlogic] = updatecurrentImage(handles,tp)

if handles.greenflag==0
    forwardlogic = 0;
    return;
else
    forwardlogic = 1;
    cellpath = handles.cellpath;
    sisterList = handles.sisterList;
    bg=handles.bg;
    row = str2num(get(handles.edit_row,'String'));
    col = str2num(get(handles.edit_col,'String'));
    field = str2num(get(handles.edit_field,'String'));
    plane = str2num(get(handles.edit_plane,'String'));
    nucCH= handles.nucCH;
    CH1= handles.CH1;
    CH2= handles.CH2;
    CH3 = handles.CH3;
    cellsize = str2num(get(handles.edit_cellsize,'String'));
    channelnames = handles.channelnames;
    totalCHs = str2num(get(handles.edit_totalCHs,'String'));
    fileformat = get(handles.edit_signalformat,'String');
    
    
    switch get(handles.popupmenu_nomin,'Value')
        case 1
            nominCH=CH1;
        case 2
            nominCH=CH2;
        case 3
            nominCH=CH3;
    end
    
    switch get(handles.popupmenu_denomin,'Value')
        case 1
            denominCH=CH1;
        case 2
            denominCH=CH2;
        case 3
            denominCH=CH3;
        case 4
            denominCH=-1;
    end
    
    if handles.filetype == 3
        filename = sprintf(fileformat,channelnames{nucCH},str2num(get(handles.edit_firstframe,'String')));
        first_info = imfinfo(fullfile(handles.ndpathname,filename));
        filename = sprintf(fileformat,channelnames{nucCH},tp);
        current_info = imfinfo(fullfile(handles.ndpathname,filename));
        [~, ~, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
        hour = 24*D+round(H);
        minute = round(MN);
        second = round(S);
        set(handles.edit_commu,'String',[num2str(hour,'%02.0f') ':' num2str(minute,'%02.0f') ':' num2str(second,'%02.0f')]);
    end
    
    filetype = handles.filetype;
    channelnames = handles.channelnames;
    template = loadsignal(handles,nucCH,tp);
    axes(handles.axes1);
    set(handles.edit_currentFrame,'String',num2str(tp));
    
    switch handles.ImageIndex
        case 3
            displayIM = loadsignal(handles,CH1,tp);
            imshow(displayIM,[str2num(get(handles.edit_ch1L,'String')) str2num(get(handles.edit_ch1H,'String'))],'Parent',handles.axes1);colormap gray;
        case 4
            displayIM = loadsignal(handles,CH2,tp);
            imshow(displayIM,[str2num(get(handles.edit_ch2L,'String')) str2num(get(handles.edit_ch2H,'String'))],'Parent',handles.axes1);colormap gray;
        case 5
            displayIM = loadsignal(handles,CH3,tp);
            imshow(displayIM,[str2num(get(handles.edit_ch3L,'String')) str2num(get(handles.edit_ch3H,'String'))],'Parent',handles.axes1);colormap gray;
        case 2
            displayIM = template;
            imshow(displayIM,[str2num(get(handles.edit_templateL,'String')) str2num(get(handles.edit_templateH,'String'))],'Parent',handles.axes1);colormap gray;
        case 1
            nominIM = loadsignal(handles,nominCH,tp);
            if denominCH == -1
                denomIM = double(ones(size(nominIM)));
            else
                denomIM = loadsignal(handles,denominCH,tp);
            end
            displayIM = calculateFRET(handles,tp,nominCH,denominCH,handles.bg);
            displayratioIm(displayIM,handles);
    end
    
    if handles.ImageIndex~=handles.overlayIndex
        
        switch handles.overlayIndex
            case 3
                displayIM = loadsignal(handles,CH1,tp);
            case 4
                displayIM = loadsignal(handles,CH2,tp);
            case 5
                displayIM = loadsignal(handles,CH3,tp);
            case 2
                displayIM = template;
            case 1
                displayIM = calculateFRET(handles,tp,nominCH,denominCH,handles.bg);
        end
        
        hold on;
        green = cat(3, zeros(size(displayIM)), ones(size(displayIM)), zeros(size(displayIM)));
        himage1 = imshow(green,'Parent',handles.axes1);hold off;
        alpha(himage1,'clear');
        
        
        
        switch handles.overlayIndex
            case 3
                alpha(himage1,mat2gray(displayIM,[str2num(get(handles.edit_ch1L,'String')) str2num(get(handles.edit_ch1H,'String'))]));
            case 4
                alpha(himage1,mat2gray(displayIM,[str2num(get(handles.edit_ch2L,'String')) str2num(get(handles.edit_ch2H,'String'))]));
            case 5
                alpha(himage1,mat2gray(displayIM,[str2num(get(handles.edit_ch3L,'String')) str2num(get(handles.edit_ch3H,'String'))]));
            case 2
                alpha(himage1,mat2gray(displayIM,[str2num(get(handles.edit_templateL,'String')) str2num(get(handles.edit_templateH,'String'))]));
            case 1
                alpha(himage1,mat2gray(displayIM,[str2num(get(handles.edit_mathL,'String')) str2num(get(handles.edit_mathH,'String'))]));
        end
        
        hold off;
        
    end
    if get(handles.checkbox_celllabeling,'Value')
        if ~isempty(cellpath)
            set(handles.listbox_cells,'String',num2str( [(1:size(cellpath{tp},1))' cellpath{tp} sisterList{tp}] ));
            
            if get(handles.checkbox_cellmarker,'Value')
                nosisterInd = find(sisterList{tp}(:,1)==-1);
                sister1Ind = find(sisterList{tp}(:,1)~=-1);
                hold on;
                plot(cellpath{tp}(nosisterInd,1),cellpath{tp}(nosisterInd,2),get(handles.edit_mark,'String'));
                plot(cellpath{tp}(sister1Ind,1),cellpath{tp}(sister1Ind,2),'o','MarkerFaceColor',get(handles.edit_sis1,'String'),'MarkerEdgeColor',get(handles.edit_sis1,'String'));
                plot(cellpath{tp}(get(handles.listbox_cells,'Value'),1),cellpath{tp}(get(handles.listbox_cells,'Value'),2),'ro');
            end
            
            if get(handles.checkbox_cellnumber,'Value')
                text(cellpath{tp}(:,1)+5,cellpath{tp}(:,2)+5,num2str((1:size(cellpath{tp},1))'),'HorizontalAlignment','left',...
                    'VerticalAlignment','middle','color',[0 .9 .5]);
            end
            hold off;
        end
        if ~isempty(bg)
            hold on;
            plot(bg{tp}(:,1),bg{tp}(:,2),'.','MarkerFaceColor','g','MarkerEdgeColor','g');
            hold off;
        end
    end
    
    drawnow;
end


function restoreCellList(handles,tp)
set(handles.listbox_cells,'Value',1);
handles.oldchosencell = 1;
set(handles.listbox_bgs,'Value',1);

cellpath = handles.cellpath;
sisterList = handles.sisterList;
cell = 1;

set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',cellpath{tp}(cell,1),cellpath{tp}(cell,2))));
set(handles.edit_cellNo,'String',num2str(cell));
set(handles.edit_sister,'String',num2str(sisterList{tp}(cell,:)));



% --- Executes on button press in pushbutton_play.
function pushbutton_play_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end
last_tp = str2num(get(handles.edit_lastframe,'String'));
c_tp = str2num(get(handles.edit_currentFrame,'String'));
cellpath = handles.cellpath;
if ~isempty(cellpath)
    last_tp = length(cellpath);
    if last_tp>str2num(get(handles.edit_lastframe,'String'));
        last_tp = str2num(get(handles.edit_lastframe,'String'));
    end
else
    last_tp = str2num(get(handles.edit_lastframe,'String'));
end

if c_tp<last_tp
    handles.greenflag = 1;
    guidata(hObject, handles);
    for tp=c_tp+1:last_tp
        handles = guidata(hObject);
        [forwardLog] = updatecurrentImage(handles,tp);
        if forwardLog==0
            break;
        end
    end
    handles.greenflag = 1;
end
guidata(hObject, handles);

function ratioIm = calculateFRET(handles,tp,nominCH,denominCH,bg)

filterParam1 = str2num(get(handles.edit_param1,'String'));
filterParam2 = str2num(get(handles.edit_param2,'String'));
bgsize = round(str2num(get(handles.edit_BGsize,'String'))/2);

signalShiftN = str2num(get(handles.edit_bg_nomin_custom,'String'));
signalShiftD = str2num(get(handles.edit_bg_denomin_custom,'String'));

nominIM = loadsignal(handles,nominCH,tp);
if denominCH == -1
    denomIM = im2double(ones(size(nominIM)));
else
    denomIM = loadsignal(handles,denominCH,tp);
end

% Load blank images if user chose to use BG points from BLANK images
if get(handles.popupmenu_bgnomin,'Value')==3
    nominBLK = loadblank(handles,nominCH,tp);
    
end
if get(handles.popupmenu_bgdenomin,'Value')==3
    if denominCH == -1
        denomBLK = im2double(ones(size(nominIM)));
    else
        denomBLK = loadblank(handles,denominCH,tp);
    end
end


if get(handles.checkbox_illumlogic,'Value')
    normN = ifft2(ifftshift(fftshift(fft2(nominIM)).*hbutter(nominIM,filterParam1,filterParam2)));
    normD = ifft2(ifftshift(fftshift(fft2(denomIM)).*hbutter(denomIM,filterParam1,filterParam2)));
else
    normN = nominIM;
    normD = denomIM;
end
BG_N = 0;
BG_D = 0;

switch get(handles.popupmenu_bgnomin,'Value')
    case 2
        if ~isempty(bg)
            for b=1:size(bg{tp},1)
                xL=max(bg{tp}(b,1)-bgsize,1);
                xR=min(bg{tp}(b,1)+bgsize,size(nominIM,2));
                yL=max(bg{tp}(b,2)-bgsize,1);
                yR=min(bg{tp}(b,2)+bgsize,size(nominIM,1));
                selectedN = normN(yL:yR,xL:xR);
                BG_N(b) = mean(selectedN(:));
            end
        end
    case 3
        xL=round(size(nominIM,2)/2)-2*bgsize;
        xR=round(size(nominIM,2)/2)+2*bgsize;
        yL=round(size(nominIM,1)/2)-2*bgsize;
        yR=round(size(nominIM,1)/2)+2*bgsize;
        BG_N = mean(nominBLK(:));
end

switch get(handles.popupmenu_bgdenomin,'Value')
    case 2
        if ~isempty(bg)
            for b=1:size(bg{tp},1)
                xL=max(bg{tp}(b,1)-bgsize,1);
                xR=min(bg{tp}(b,1)+bgsize,size(nominIM,2));
                yL=max(bg{tp}(b,2)-bgsize,1);
                yR=min(bg{tp}(b,2)+bgsize,size(nominIM,1));
                selectedD = normD(yL:yR,xL:xR);
                BG_D(b) = mean(selectedD(:));
            end
        end
    case 3
        xL=round(size(nominIM,2)/2)-2*bgsize;
        xR=round(size(nominIM,2)/2)+2*bgsize;
        yL=round(size(nominIM,1)/2)-2*bgsize;
        yR=round(size(nominIM,1)/2)+2*bgsize;
        BG_D = mean(denomBLK(:));
end


normN = normN-mean(BG_N);
normD = normD-mean(BG_D);



normN = normN+signalShiftN;
normD = normD+signalShiftD;

ratioIm =  normN./normD;


function displayratioIm(ratioIm,handles)
displaygateL = str2num(get(handles.edit_mathL,'String'));
displaygateH = str2num(get(handles.edit_mathH,'String'));

imshow(ratioIm,[displaygateL displaygateH],'Parent',handles.axes1);hold on;
switch get(handles.popupmenu_colormap,'Value')
    case 1
        set(gcf,'Colormap',handles.mycmap1)
    case 2
        set(gcf,'Colormap',handles.mycmap2)
    case 3
        set(gcf,'Colormap',handles.mycmap3)
    case 4
        set(gcf,'Colormap',handles.mycmap4)
end



% --- Executes on button press in pushbutton_pause.
function pushbutton_pause_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.greenflag = 0;
guidata(hObject, handles);

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
        set(handles.edit_commu,'String','Cell edging (Canny) chosen');
        handles.cellrecogmethod = 1;
    case 'radiobutton_method2'
        set(handles.edit_commu,'String','Contour tracking chosen');
        handles.cellrecogmethod = 2;
    otherwise
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles) %#ok<*INUSL>
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

load(fullfile(handles.celltrackpathname,get(handles.edit_outputname,'String')));

tp=str2num(get(handles.edit_currentFrame,'String'));

handles.cellpath=cellpath;
handles.sisterList=sisterList;
if exist('bg','var')
    handles.bg = bg;
else
    handles.bg = [];
end
guidata(hObject, handles);
handles = guidata(hObject);
updatecurrentImage(handles,tp);


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mysignal = handles.mysignal;
savefile = ['AllFRETsignals_r' num2str(row) 'c' num2str(col) 'f' num2str(field)];
tpframe = ( str2num(get(handles.edit_minstamp,'String')) + str2num(get(handles.edit_secstamp,'String'))/60 ); %minutes

first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

timestep = tpframe*((first_tp:last_tp)-1); %minutes

save(savefile,'mysignal','timestep','-v7.3')


% --- Executes during object creation, after setting all properties.
function listbox_cells_CreateFcn(hObject, ~, ~)
% hObject    handle to listbox_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

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

function edit_nucCH_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nucCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nucCH as text
%        str2double(get(hObject,'String')) returns contents of edit_nucCH as a double
handles.nucCH = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_nucCH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nucCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.nucCH= str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_plane_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
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
%        str2double(get(hObject,'String')) returns contents of edit_lastframe as a double


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
function edit_sigma_CreateFcn(hObject, ~, handles)
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


function edit_maxI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxI as text
%        str2double(get(hObject,'String')) returns contents of edit_maxI as a double
handles.maxI = str2num(get(hObject,'String'));
guidata(hObject, handles);

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



function edit_CH2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CH2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CH2 as text
%        str2double(get(hObject,'String')) returns contents of edit_CH2 as a double


% --- Executes during object creation, after setting all properties.
function edit_CH2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CH2 (see GCBO)
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



function edit_signalformat_Callback(hObject, eventdata, handles)
% hObject    handle to edit_signalformat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_signalformat as text
%        str2double(get(hObject,'String')) returns contents of edit_signalformat as a double


% --- Executes during object creation, after setting all properties.
function edit_signalformat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_signalformat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_label.
function checkbox_label_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_label


% --- Executes on button press in pushbutton_genim.
function pushbutton_genim_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_genim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
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
nucCH= str2num(get(handles.edit_nucCH,'String'));
CH1= handles.CH1;
CH2= handles.CH2;
CH3 = handles.CH3;
filterParam = [2 2];
cellsize = str2num(get(handles.edit_cellsize,'String'));
totalCHs = str2num(get(handles.edit_totalCHs,'String'));
fileformat = get(handles.edit_signalformat,'String');

displaygateL = str2num(get(handles.edit_mathL,'String'));
displaygateH = str2num(get(handles.edit_mathH,'String'));
tp=str2num(get(handles.edit_currentFrame,'String'));

first_tp = str2num(get(handles.edit_firstframe,'String'));
if handles.filetype == 3
    filename = sprintf(fileformat,handles.channelnames{nucCH},first_tp);
    first_info = imfinfo(fullfile(handles.ndpathname,filename));
end


switch get(handles.popupmenu_nomin,'Value')
    case 1
        nominCH=CH1;
    case 2
        nominCH=CH2;
    case 3
        nominCH=CH3;
end

switch get(handles.popupmenu_denomin,'Value')
    case 1
        denominCH=CH1;
    case 2
        denominCH=CH2;
    case 3
        denominCH=CH3;
    case 4
        denominCH=-1;
end
figure();
outer1 = gca;


filetype = handles.filetype;
channelnames = handles.channelnames;
template = loadsignal(handles,nucCH,tp);
if handles.filetype == 3
    filename = sprintf(fileformat,handles.channelnames{nucCH},tp);
    current_info = imfinfo(fullfile(handles.ndpathname,filename));
end
switch handles.ImageIndex
    case 3
        displayIM = loadsignal(handles,CH1,tp);
        imshow(displayIM,[str2num(get(handles.edit_ch1L,'String')) str2num(get(handles.edit_ch1H,'String'))],'Parent',gca);colormap(gray);
    case 4
        displayIM = loadsignal(handles,CH2,tp);
        imshow(displayIM,[str2num(get(handles.edit_ch2L,'String')) str2num(get(handles.edit_ch2H,'String'))],'Parent',gca);colormap(gray);
    case 5
        displayIM = loadsignal(handles,CH3,tp);
        imshow(displayIM,[str2num(get(handles.edit_ch3L,'String')) str2num(get(handles.edit_ch3H,'String'))],'Parent',gca);colormap(gray);
    case 2
        displayIM = template;
        imshow(displayIM,[str2num(get(handles.edit_templateL,'String')) str2num(get(handles.edit_templateH,'String'))],'Parent',gca);colormap(gray);
    case 1
        nominIM = loadsignal(handles,nominCH,tp);
        displayIM = calculateFRET(handles,tp,nominCH,denominCH,[]);
        imshow(displayIM,[str2num(get(handles.edit_mathL,'String')) str2num(get(handles.edit_mathH,'String'))],'Parent',gca);
        switch get(handles.popupmenu_colormap,'Value')
            case 1
                set(gcf,'Colormap',handles.mycmap1)
            case 2
                set(gcf,'Colormap',handles.mycmap2)
            case 3
                set(gcf,'Colormap',handles.mycmap3)
            case 4
                set(gcf,'Colormap',handles.mycmap4)
        end
end

hold on;
if ~isempty(cellpath)
    maxCell = size(cellpath{end},1);
    rancolor = lines(maxCell);
    for cell=1:size(cellpath{tp},1)
        if sisterList{tp}(cell)==-1
            plot(cellpath{tp}(cell,1),cellpath{tp}(cell,2),'rx');
        else
            plot(cellpath{tp}(cell,1),cellpath{tp}(cell,2),'o','MarkerFaceColor',rancolor(1+mod(cell,maxCell),:),'MarkerEdgeColor',rancolor(1+mod(cell,maxCell),:));
            plot(cellpath{tp}(sisterList{tp}(cell),1),cellpath{tp}(sisterList{tp}(cell),2),'o','MarkerFaceColor',rancolor(1+mod(cell,maxCell),:),'MarkerEdgeColor',rancolor(1+mod(cell,maxCell),:));
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
                [~, ~, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
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


% --- Executes on button press in pushbutton_genmov.
function pushbutton_genmov_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_genmov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
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
nucCH= str2num(get(handles.edit_nucCH,'String'));
CH1= handles.CH1;
CH2= handles.CH2;
CH3 = handles.CH3;
filterParam = [2 2];
cellsize = str2num(get(handles.edit_cellsize,'String'));
totalCHs = str2num(get(handles.edit_totalCHs,'String'));
fileformat = get(handles.edit_signalformat,'String');

displaygateL = str2num(get(handles.edit_mathL,'String'));
displaygateH = str2num(get(handles.edit_mathH,'String'));

first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));

switch get(handles.popupmenu_nomin,'Value')
    case 1
        nominCH=CH1;
    case 2
        nominCH=CH2;
    case 3
        nominCH=CH3;
end

switch get(handles.popupmenu_denomin,'Value')
    case 1
        denominCH=CH1;
    case 2
        denominCH=CH2;
    case 3
        denominCH=CH3;
    case 4
        denominCH=-1;
end

if ~isempty(bg)
    allbg = cat(1,bg{:});
    imshift(:,1) =round( -mean(allbg(:,1)) + allbg(:,1) );
    imshift(:,2) =round( -mean(allbg(:,2)) + allbg(:,2) );
    maxshiftX = max(abs(imshift(:,1)));
    maxshiftY = max(abs(imshift(:,2)));
end

aviobj = VideoWriter(fullfile(handles.ndpathname,['myMov_r' num2str(row) 'c' num2str(col) 'f' num2str(field) 'ch' num2str(handles.ImageIndex) '.avi']));
aviobj.FrameRate = str2num(get(handles.edit_framerate,'String'));
aviobj.Quality = 80;
open(aviobj);
f1=figure();
axes(gca);
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

if handles.filetype == 3
    filename = sprintf(fileformat,handles.channelnames{nucCH},first_tp);
    first_info = imfinfo(fullfile(handles.ndpathname,filename));
end


for tp=first_tp:last_tp
    
    
    if handles.ImageIndex == handles.overlayIndex
        
        template = loadsignal(handles,nucCH,tp);
        if handles.filetype == 3
            filename = sprintf(fileformat,handles.channelnames{nucCH},tp);
            current_info = imfinfo(fullfile(handles.ndpathname,filename));
        end
        
        frameX = [1 size(template,2)];
        frameY = [1 size(template,1)];
        
        if ~isempty(bg) && get(handles.framshift_logic,'Value')
            frameX =  [ maxshiftX size(template,2)-maxshiftX ] + imshift(tp,1) ;
            frameY =  [ maxshiftY size(template,1)-maxshiftY ] + imshift(tp,2) ;
            
            if frameX(1) <= 0
                frameX(1) = 1;
                frameX(2) = frameX(2)+1;
            end
            
            if frameX(2) > size(template,2)
                frameX(2) = size(template,2);
            end
            %+ [-min(imshift(:,1))   size(template,2) - max(imshift(:,1)) ];
            
            if frameY(1) <= 0
                frameY(1) = 1;
                frameY(2) = frameY(2)+1;
            end
            
            if frameY(2) > size(template,1)
                frameY(2) = size(template,1);
            end
            %+ [-min(imshift(:,2))   size(template,1) - max(imshift(:,2)) ];
        end
        
        switch handles.ImageIndex
            
            case 3
                displayIM = loadsignal(handles,CH1,tp);
                displayIM = displayIM(frameY(1):frameY(2),frameX(1):frameX(2));
                
                imshow(displayIM,[str2num(get(handles.edit_ch1L,'String')) str2num(get(handles.edit_ch1H,'String'))],'Parent',gca,'InitialMagnification',50);
                colormap(gray);
            case 4
                displayIM = loadsignal(handles,CH2,tp);
                displayIM = displayIM(frameY(1):frameY(2),frameX(1):frameX(2));
                
                imshow(displayIM,[str2num(get(handles.edit_ch2L,'String')) str2num(get(handles.edit_ch2H,'String'))],'Parent',gca,'InitialMagnification',50);
                colormap(gray);
            case 5
                displayIM = loadsignal(handles,CH3,tp);
                displayIM = displayIM(frameY(1):frameY(2),frameX(1):frameX(2));
                
                imshow(displayIM,[str2num(get(handles.edit_ch3L,'String')) str2num(get(handles.edit_ch3H,'String'))],'Parent',gca,'InitialMagnification',50);
                colormap(gray);
            case 2
                displayIM = template;
                displayIM = displayIM(frameY(1):frameY(2),frameX(1):frameX(2));
                
                imshow(displayIM,[str2num(get(handles.edit_templateL,'String')) str2num(get(handles.edit_templateH,'String'))],'Parent',gca,'InitialMagnification',50);
                colormap(gray);
            case 1
                displayIM = calculateFRET(handles,tp,nominCH,denominCH,bg);
                displayIM = displayIM(frameY(1):frameY(2),frameX(1):frameX(2));
                
                imshow(displayIM,[str2num(get(handles.edit_mathL,'String')) str2num(get(handles.edit_mathH,'String'))],'Parent',gca,'InitialMagnification',50);
                switch get(handles.popupmenu_colormap,'Value')
                    case 1
                        set(gcf,'Colormap',handles.mycmap1)
                    case 2
                        set(gcf,'Colormap',handles.mycmap2)
                    case 3
                        set(gcf,'Colormap',handles.mycmap3)
                    case 4
                        set(gcf,'Colormap',handles.mycmap4)
                end
        end
        
    else
        
        template = loadsignal(handles,nucCH,tp);
        if handles.filetype == 3
            filename = sprintf(fileformat,handles.channelnames{nucCH},tp);
            current_info = imfinfo(fullfile(handles.ndpathname,filename));
        end
        switch handles.ImageIndex
            case 3
                I1 = mat2gray(loadsignal(handles,CH1,tp),[str2num(get(handles.edit_ch1L,'String')) str2num(get(handles.edit_ch1H,'String'))]);
            case 4
                I1 = mat2gray(loadsignal(handles,CH2,tp),[str2num(get(handles.edit_ch2L,'String')) str2num(get(handles.edit_ch2H,'String'))]);
            case 5
                I1 = mat2gray(loadsignal(handles,CH3,tp),[str2num(get(handles.edit_ch3L,'String')) str2num(get(handles.edit_ch3H,'String'))]);
            case 2
                I1 = mat2gray(template,[str2num(get(handles.edit_templateL,'String')) str2num(get(handles.edit_templateH,'String'))]);
            case 1
                I1 = mat2gray(calculateFRET(handles,tp,nominCH,denominCH,bg),[str2num(get(handles.edit_mathL,'String')) str2num(get(handles.edit_mathH,'String'))]);
        end
        
        switch handles.overlayIndex
            
            case 3
                I2 = mat2gray(loadsignal(handles,CH1,tp),[str2num(get(handles.edit_ch1L,'String')) str2num(get(handles.edit_ch1H,'String'))]);
                
            case 4
                I2 = mat2gray(loadsignal(handles,CH2,tp),[str2num(get(handles.edit_ch2L,'String')) str2num(get(handles.edit_ch2H,'String'))]);
                
            case 5
                I2 = mat2gray(loadsignal(handles,CH3,tp),[str2num(get(handles.edit_ch3L,'String')) str2num(get(handles.edit_ch3H,'String'))]);
                
            case 2
                I2 = mat2gray(template,[str2num(get(handles.edit_templateL,'String')) str2num(get(handles.edit_templateH,'String'))]);
                
            case 1
                I2 = mat2gray(calculateFRET(handles,tp,nominCH,denominCH,bg),[str2num(get(handles.edit_mathL,'String')) str2num(get(handles.edit_mathH,'String'))]);
        end
        
        combinedI(:,:,1) = 0.8*I1;
        combinedI(:,:,2) = 0.8*I1+0.8*I2;
        combinedI(:,:,3) = 0.8*I1;
        imshow(combinedI,'Parent',gca,'InitialMagnification',50);
        
    end
    
    hold on;
    
    if ~isempty(cellpath) && ~isempty(bg)
        
        
        maxCell = size(cellpath{end},1);
        rancolor = lines(maxCell);
        if get(handles.framshift_logic,'Value')
            
            for cell=1:size(cellpath{tp},1)
                if(cellpath{tp}(cell,1)-frameX(1) > frameX(1)+20 && cellpath{tp}(cell,1)-frameX(1) < frameX(2)-20 && ...
                        cellpath{tp}(cell,2)-frameY(1) > frameY(1)+20 && cellpath{tp}(cell,2)-frameY(1) < frameY(2)-20)
                    if sisterList{tp}(cell)==-1
                        plot(cellpath{tp}(cell,1)-frameX(1),cellpath{tp}(cell,2)-frameY(1),'rx');
                    else
                        plot(cellpath{tp}(cell,1)-frameX(1),cellpath{tp}(cell,2)-frameY(1),'o','MarkerFaceColor',rancolor(1+mod(cell,maxCell),:),'MarkerEdgeColor',rancolor(1+mod(cell,maxCell),:));
                        plot(cellpath{tp}(sisterList{tp}(cell),1)-frameX(1),cellpath{tp}(sisterList{tp}(cell),2)-frameY(1),'o','MarkerFaceColor',rancolor(1+mod(cell,maxCell),:),'MarkerEdgeColor',rancolor(1+mod(cell,maxCell),:));
                    end
                    
                    text(cellpath{tp}(cell,1)+10-frameX(1),cellpath{tp}(cell,2)-10-frameY(1),num2str(cell),'HorizontalAlignment','left',...
                        'VerticalAlignment','middle','color',[0 .9 .5]);
                end
            end
            plot(bg{tp}(:,1)-frameX(1), bg{tp}(:,2)-frameY(1),'ko','MarkerFaceColor','k');
        else
            for cell=1:size(cellpath{tp},1)
                if sisterList{tp}(cell)==-1
                    
                    plot(cellpath{tp}(cell,1),cellpath{tp}(cell,2),'rx');
                else
                    plot(cellpath{tp}(cell,1),cellpath{tp}(cell,2),'o','MarkerFaceColor',rancolor(1+mod(cell,maxCell),:),'MarkerEdgeColor',rancolor(1+mod(cell,maxCell),:));
                    plot(cellpath{tp}(sisterList{tp}(cell),1),cellpath{tp}(sisterList{tp}(cell),2),'o','MarkerFaceColor',rancolor(1+mod(cell,maxCell),:),'MarkerEdgeColor',rancolor(1+mod(cell,maxCell),:));
                end
                
                text(cellpath{tp}(cell,1)+10,cellpath{tp}(cell,2)+10,num2str(cell),'HorizontalAlignment','left',...
                    'VerticalAlignment','middle','color',[0 .9 .5]);
            end
            
        end
        
    end
    switch handles.framestamp
        case 1
            text(30,30,['frame:' num2str(tp)],...
                'HorizontalAlignment','left','VerticalAlignment','top','color','w','BackgroundColor','k','fontsize',22);
        case 2
            switch handles.filetype
                case 3
                    [~, ~, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
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
    axis image;
    drawnow;
    F = getframe(gca);
    writeVideo(aviobj,F);
    
end
close(aviobj);
close(f1);

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
        handles.filetype = 1;
    case 'radiobutton_tiffstack'
        handles.filetype = 2;
    case 'radiobutton_customtiff'
        handles.filetype = 3;
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton_locatefile.
function pushbutton_locatefile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_locatefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname,FilterIndex] = uigetfile(...
    {'*.tif','stacked tiff Files (*.tif)';
    '*.*',  'All Files (*.*)'}, 'Pick a file');
if FilterIndex~=0
    set(handles.edit_signalformat,'String',filename);
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
handles.filetype = 1;
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


function edit_fret_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fret (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fret as text
%        str2double(get(hObject,'String')) returns contents of edit_fret as a double

handles.CH1= str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_fret_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fret (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.CH1= str2num(get(hObject,'String'));
guidata(hObject, handles);


function edit_cfp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cfp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cfp as text
%        str2double(get(hObject,'String')) returns contents of edit_cfp as a double
handles.CH2= str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_cfp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cfp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.CH2= str2num(get(hObject,'String'));
guidata(hObject, handles);

function edit_yfp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_yfp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_yfp as text
%        str2double(get(hObject,'String')) returns contents of edit_yfp as a double

handles.CH3= str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_yfp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_yfp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.CH3= str2num(get(hObject,'String'));
guidata(hObject, handles);

function edit_outputname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_outputname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_outputname as text
%        str2double(get(hObject,'String')) returns contents of edit_outputname as a double


% --- Executes during object creation, after setting all properties.
function edit_outputname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_outputname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_celltracklocator.

% --- Executes on selection change in popupmenu_nomin.
function popupmenu_nomin_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_nomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_nomin contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_nomin


% --- Executes during object creation, after setting all properties.
function popupmenu_nomin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_nomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_bgnomin.
function popupmenu_bgnomin_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_bgnomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_bgnomin contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_bgnomin

switch get(hObject,'Value')
    case 1
        set(handles.edit_bg_nomin_custom,'String','0');
    case 2
        set(handles.edit_bg_nomin_custom,'String','0.001');
    case 3
        set(handles.edit_bg_nomin_custom,'String','0.001');
end
% --- Executes during object creation, after setting all properties.
function popupmenu_bgnomin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_bgnomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_denomin.
function popupmenu_denomin_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_denomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_denomin contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_denomin


% --- Executes during object creation, after setting all properties.
function popupmenu_denomin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_denomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_bgdenomin.
function popupmenu_bgdenomin_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_bgdenomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_bgdenomin contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_bgdenomin
switch get(hObject,'Value')
    case 1
        set(handles.edit_bg_denomin_custom,'String','0');
    case 2
        set(handles.edit_bg_denomin_custom,'String','0.001');
    case 3
        set(handles.edit_bg_denomin_custom,'String','0.001');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_bgdenomin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_bgdenomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_illumlogic.
function checkbox_illumlogic_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_illumlogic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_illumlogic
switch get(hObject,'Value')
    case get(hObject,'Max')
        set(handles.edit_bg_nomin_custom,'String','0.001');
        set(handles.edit_bg_denomin_custom,'String','0.001');
    case get(hObject,'Min')
        set(handles.edit_bg_nomin_custom,'String','0');
        set(handles.edit_bg_denomin_custom,'String','0');
end
function edit_mathL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mathL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mathL as text
%        str2double(get(hObject,'String')) returns contents of edit_mathL as a double


% --- Executes during object creation, after setting all properties.
function edit_mathL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mathL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mathH_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mathH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mathH as text
%        str2double(get(hObject,'String')) returns contents of edit_mathH as a double

% --- Executes during object creation, after setting all properties.
function edit_mathH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mathH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bg_nomin_custom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bg_nomin_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bg_nomin_custom as text
%        str2double(get(hObject,'String')) returns contents of edit_bg_nomin_custom as a double


% --- Executes during object creation, after setting all properties.
function edit_bg_nomin_custom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bg_nomin_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bg_denomin_custom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bg_denomin_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bg_denomin_custom as text
%        str2double(get(hObject,'String')) returns contents of edit_bg_denomin_custom as a double


% --- Executes during object creation, after setting all properties.
function edit_bg_denomin_custom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bg_denomin_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_colormap.
function popupmenu_colormap_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_colormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_colormap
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

tp=str2num(get(handles.edit_currentFrame,'String'));
updatecurrentImage(handles,tp);
guidata(hObject, handles);

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


% --- Executes on button press in pushbutton_collect.
function pushbutton_collect_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to pushbutton_collect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CH1= handles.CH1;
CH2= handles.CH2;
CH3 = handles.CH3;
switch get(handles.popupmenu_nomin,'Value')
    case 1
        nominCH=CH1;
    case 2
        nominCH=CH2;
    case 3
        nominCH=CH3;
end
handles.nominCH= nominCH;

switch get(handles.popupmenu_denomin,'Value')
    case 1
        denominCH=CH1;
    case 2
        denominCH=CH2;
    case 3
        denominCH=CH3;
    case 4
        denominCH=-1;
end
handles.denominCH= denominCH;
handles.nucCH = str2num(get(handles.edit_nucCH,'String'));
handles.cellCH = str2num(get(handles.edit_cellCH,'String'));
handles.cellsize = str2num(get(handles.edit_cellsize,'String'));
handles.totalCHs = str2num(get(handles.edit_totalCHs,'String'));
handles.signalformat = get(handles.edit_signalformat,'String');
handles.blankformat = get(handles.edit_fileformatBG,'String');
handles.first_tp = str2num(get(handles.edit_firstframe,'String'));
handles.last_tp = str2num(get(handles.edit_lastframe,'String'));
handles.minstamp = get(handles.edit_minstamp,'String');
handles.secstamp = get(handles.edit_secstamp,'String');
handles.cytosize = get(handles.edit_cytosize,'String');

handles.bgnominType = get(handles.popupmenu_bgnomin,'Value');
handles.bgdenominType = get(handles.popupmenu_bgdenomin,'Value');
handles.illumlogic = get(handles.checkbox_illumlogic,'Value');

handles.Var1LOG = get(handles.checkbox_variable1,'Value');
handles.Var2LOG = get(handles.checkbox_variable2,'Value');
handles.Var3LOG = get(handles.checkbox_variable3,'Value');
handles.Var4LOG = get(handles.checkbox_variable4,'Value');

handles.region1LOG = get(handles.popupmenu_regionVar1,'Value');
handles.region2LOG = get(handles.popupmenu_regionVar2,'Value');
handles.region3LOG = get(handles.popupmenu_regionVar3,'Value');
handles.region4LOG = get(handles.popupmenu_regionVar4,'Value');

handles.signal1LOG = get(handles.popupmenu_signal1,'Value');
handles.signal2LOG = get(handles.popupmenu_signal2,'Value');
handles.signal3LOG = get(handles.popupmenu_signal3,'Value');
handles.signal4LOG = get(handles.popupmenu_signal4,'Value');

handles.regions1 = get(handles.popupmenu_regionVar1,'String');
handles.regions2 = get(handles.popupmenu_regionVar2,'String');
handles.regions3 = get(handles.popupmenu_regionVar3,'String');
handles.regions4 = get(handles.popupmenu_regionVar4,'String');

handles.signals1 = get(handles.popupmenu_signal1,'String');
handles.signals2 = get(handles.popupmenu_signal2,'String');
handles.signals3 = get(handles.popupmenu_signal3,'String');
handles.signals4 = get(handles.popupmenu_signal4,'String');

handles.outputsignaldatasetname = get(handles.edit_outputname,'String');

handles.filterParam1 = str2num(get(handles.edit_param1,'String'));
handles.filterParam2 = str2num(get(handles.edit_param2,'String'));
handles.bgsize = round(str2num(get(handles.edit_BGsize,'String'))/2);

handles.signalShiftN = str2num(get(handles.edit_bg_nomin_custom,'String'));
handles.signalShiftD = str2num(get(handles.edit_bg_denomin_custom,'String'));

guidata(hObject, handles);
handles = guidata(hObject);

signalformat = get(handles.edit_signalformat,'String');
blankformat = get(handles.edit_fileformatBG,'String');
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
collectdata_individual(handles,row,col,field,signalformat,blankformat);
updateOutputList(handles);

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
    
    SisList = unique(thirdSis);
    
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


function updateOutputList(handles)



row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
if exist(fullfile(handles.ndpathname,H5filename),'file')
    
    info = h5info(fullfile(handles.ndpathname,H5filename),['/field' num2str(field)]);
    basename='signal';
    cindex = 1;
    outputcounts=[];
    outputList=[];
    if length(info.Datasets)>0 %#ok<*ISMT>
        for i=1:length(info.Datasets)
            tmp = regexp(info.Datasets(i).Name,basename, 'match');
            if ~isempty(tmp)
                outputList{cindex}= info.Datasets(i).Name;
                cindex=cindex+1;
            end
        end
    end
    if isempty(outputList)
        outputList = '<not present>';
    end
    set(handles.popupmenu_output,'String',outputList);
    set(handles.popupmenu_output,'Value',1);
end
function outputsignal = signalOutputing(regiontype,signaltype,mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y)
nucData = [];
cytoData = [];
cellData = [];
switch regiontype
    case 1
        switch signaltype
            case 1
                for i=1:length(nuc_Y)
                    nucData(i) = mini_ratioIm(nuc_Y(i),nuc_X(i));
                end
                for i=1:length(cyto_Y)
                    cytoData(i) = mini_ratioIm(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(nucData)
                    Average_Nuc = mean(nucData);
                else
                    Average_Nuc = 0;
                end
                if ~isempty(cytoData)
                    Average_Cyto = mean(cytoData);
                else
                    Average_Cyto = 0;
                end
            case 2
                for i=1:length(nuc_Y)
                    nucData(i) = cell_ch1(nuc_Y(i),nuc_X(i));
                end

                for i=1:length(cyto_Y)
                    cytoData(i) = cell_ch1(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(nucData)
                    Average_Nuc = mean(nucData);
                else
                    Average_Nuc = 0;
                end
                if ~isempty(cytoData)
                    Average_Cyto = mean(cytoData);
                else
                    Average_Cyto = 0;
                end
                
            case 3
                for i=1:length(nuc_Y)
                    nucData(i) = cell_ch2(nuc_Y(i),nuc_X(i));
                end
                for i=1:length(cyto_Y)
                    cytoData(i) = cell_ch2(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(nucData)
                    Average_Nuc = mean(nucData);
                else
                    Average_Nuc = 0;
                end
                if ~isempty(cytoData)
                    Average_Cyto = mean(cytoData);
                else
                    Average_Cyto = 0;
                end
            case 4
                for i=1:length(nuc_Y)
                    nucData(i) = cell_ch3(nuc_Y(i),nuc_X(i));
                end

                for i=1:length(cyto_Y)
                    cytoData(i) = cell_ch3(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(nucData)
                    Average_Nuc = mean(nucData);
                else
                    Average_Nuc = 0;
                end
                if ~isempty(cytoData)
                    Average_Cyto = mean(cytoData);
                else
                    Average_Cyto = 0;
                end
        end
        if Average_Cyto == 0 || Average_Nuc == 0
            outputsignal = 0;
        else
            outputsignal = Average_Nuc/Average_Cyto;
        end
    case 2
        switch signaltype
            case 1
                for i=1:length(nuc_Y)
                    nucData(i) = mini_ratioIm(nuc_Y(i),nuc_X(i));
                end
                if ~isempty(nucData)
                    outputsignal = mean(nucData);
                else
                    outputsignal = 0;
                end
            case 2
                for i=1:length(nuc_Y)
                    nucData(i) = cell_ch1(nuc_Y(i),nuc_X(i));
                end
                if ~isempty(nucData)
                    outputsignal = mean(nucData);
                else
                    outputsignal = 0;
                end
            case 3
                for i=1:length(nuc_Y)
                    nucData(i) = cell_ch2(nuc_Y(i),nuc_X(i));
                end
                if ~isempty(nucData)
                    outputsignal = mean(nucData);
                else
                    outputsignal = 0;
                end
            case 4
                for i=1:length(nuc_Y)
                    nucData(i) = cell_ch3(nuc_Y(i),nuc_X(i));
                end
                if ~isempty(nucData)
                    outputsignal = mean(nucData);
                else
                    outputsignal = 0;
                end
        end
    case 3
        switch signaltype
            case 1
                for i=1:length(cyto_Y)
                    cytoData(i) = mini_ratioIm(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(cytoData)
                    outputsignal = mean(cytoData);
                else
                    outputsignal = 0;
                end
            case 2
                for i=1:length(cyto_Y)
                    cytoData(i) = cell_ch1(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(cytoData)
                    outputsignal = mean(cytoData);
                else
                    outputsignal = 0;
                end
            case 3
                for i=1:length(cyto_Y)
                    cytoData(i) = cell_ch2(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(cytoData)
                    outputsignal = mean(cytoData);
                else
                    outputsignal = 0;
                end
            case 4
                for i=1:length(cyto_Y)
                    cytoData(i) = cell_ch3(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(cytoData)
                    outputsignal = mean(cytoData);
                else
                    outputsignal = 0;
                end
        end
    case 4
        switch signaltype
            case 1
                for i=1:length(cell_Y)
                    cellData(i) = mini_ratioIm(cell_Y(i),cell_X(i));
                end
                if ~isempty(cellData)
                    outputsignal = mean(cellData);
                else
                    outputsignal = 0;
                end
            case 2
                for i=1:length(cell_Y)
                    cellData(i) = cell_ch1(cell_Y(i),cell_X(i));
                end
                if ~isempty(cellData)
                    outputsignal = mean(cellData);
                else
                    outputsignal = 0;
                end
            case 3
                for i=1:length(cell_Y)
                    cellData(i) = cell_ch2(cell_Y(i),cell_X(i));
                end
                if ~isempty(cellData)
                    outputsignal = mean(cellData);
                else
                    outputsignal = 0;
                end
            case 4
                for i=1:length(cell_Y)
                    cellData(i) = cell_ch3(cell_Y(i),cell_X(i));
                end
                if ~isempty(cellData)
                    outputsignal = mean(cellData);
                else
                    outputsignal = 0;
                end
        end
end

function edit_param1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_param1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_param1 as text
%        str2double(get(hObject,'String')) returns contents of edit_param1 as a double


% --- Executes during object creation, after setting all properties.
function edit_param1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_param1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_param2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_param2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_param2 as text
%        str2double(get(hObject,'String')) returns contents of edit_param2 as a double


% --- Executes during object creation, after setting all properties.
function edit_param2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_param2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

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

function edit_sis1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sis1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sis1 as text
%        str2double(get(hObject,'String')) returns contents of edit_sis1 as a double


% --- Executes during object creation, after setting all properties.
function edit_sis1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sis1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [BW, BWperim] = CentroidToMask2(M,xg,yg,maxI,thresLow,thresHigh,sigma,invertLog)
if invertLog
    inverted = ((maxI)-1)-M;
    BW = imfill(edge(inverted,'canny',[thresLow thresHigh],sigma),'holes');
else
    BW = imfill(edge(M,'canny',[thresLow thresHigh],sigma),'holes');
end

BW = bwselect(BW,xg,yg);
BW2 = bwmorph(BW,'dilate',1);
BW3 = bwmorph(BW,'dilate',4);
BW4 = BW3-BW2;
BWperim = bwmorph(BW,'remove','inf');



function [BW4, BWperim]= CentroidToMask(M,xg,yg,filterParam,sizeThres)

inverted = (((2^15)-1)-M);

thres = linspace(0.001,0.002,4);
oldBW = ones(size(inverted));
for i=1:4
    BW = edge(inverted,'sobel',thres(i),'nothinning');
    se = strel('disk',3);
    closeBW = imclose(BW,se);
    closeBW = bwmorph(closeBW,'remove',Inf);
    im = bwselect(imfill(closeBW,[xg yg]),xg,yg);
    mask{i} = bwmorph(im,'majority',10);
    cellarea(i) = bwarea(mask{i});
    
    
end
drawnow;
ind = find(cellarea>sizeThres(1) & cellarea<sizeThres(2));
if ~isempty(ind)
    outind = find(max(cellarea(ind))==cellarea);
    outBW = mask{outind};
    BW2 = bwmorph(outBW,'erode',1);
    BW3 = bwmorph(outBW,'dilate',3);
    BW4 = BW3-BW2;
else
    BW4 = im2bw(M,0.025);
end
BWperim = bwmorph(BW4,'remove',Inf);


% --- Executes during object creation, after setting all properties.
function pushbutton_plotindividual_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_plotindividual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_colormapeditor.
function pushbutton_colormapeditor_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_colormapeditor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

axes(handles.axes1);
colormapeditor;


% --- Executes on button press in pushbutton_savecolormap.
function pushbutton_savecolormap_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_savecolormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

mycmap = get(gcf,'Colormap');
assignin('base','mycmap',mycmap);
switch get(handles.popupmenu_colormap,'Value')
    case 1
        handles.mycmap1=mycmap;
    case 2
        handles.mycmap2=mycmap;
    case 3
        handles.mycmap3=mycmap;
    case 4
        handles.mycmap4=mycmap;
    case 5
        handles.mycmap5=mycmap;
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function uipanel_datacollectiontype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_datacollectiontype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.collectiontypeind = 1;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_cytosize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cytosize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_xshift_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xshift as text
%        str2double(get(hObject,'String')) returns contents of edit_xshift as a double


% --- Executes during object creation, after setting all properties.
function edit_xshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function edit_collectedcells_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_collectedcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function edit_framerate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_regionVar1.
function popupmenu_regionVar1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_regionVar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_regionVar1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_regionVar1


% --- Executes during object creation, after setting all properties.
function popupmenu_regionVar1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_regionVar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_signal1.
function popupmenu_signal1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_signal1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_signal1


% --- Executes during object creation, after setting all properties.
function popupmenu_signal1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_variable1.
function checkbox_variable1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_variable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_variable1


% --- Executes on selection change in popupmenu_regionVar2.
function popupmenu_regionVar2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_regionVar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_regionVar2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_regionVar2


% --- Executes during object creation, after setting all properties.
function popupmenu_regionVar2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_regionVar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_signal2.
function popupmenu_signal2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_signal2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_signal2


% --- Executes during object creation, after setting all properties.
function popupmenu_signal2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_variable2.
function checkbox_variable2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_variable2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_variable2


% --- Executes on selection change in popupmenu_regionVar3.
function popupmenu_regionVar3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_regionVar3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_regionVar3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_regionVar3


% --- Executes during object creation, after setting all properties.
function popupmenu_regionVar3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_regionVar3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_signal3.
function popupmenu_signal3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_signal3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_signal3


% --- Executes during object creation, after setting all properties.
function popupmenu_signal3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_variable3.
function checkbox_variable3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_variable3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_variable3


% --- Executes on selection change in popupmenu_regionVar4.
function popupmenu_regionVar4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_regionVar4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_regionVar4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_regionVar4


% --- Executes during object creation, after setting all properties.
function popupmenu_regionVar4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_regionVar4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_signal4.
function popupmenu_signal4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_signal4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_signal4


% --- Executes during object creation, after setting all properties.
function popupmenu_signal4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_variable4.
function checkbox_variable4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_variable4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_variable4


function edit_var1name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_var1name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_var1name as text
%        str2double(get(hObject,'String')) returns contents of edit_var1name as a double


% --- Executes during object creation, after setting all properties.
function edit_var1name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_var1name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_var2name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_var2name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_var2name as text
%        str2double(get(hObject,'String')) returns contents of edit_var2name as a double


% --- Executes during object creation, after setting all properties.
function edit_var2name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_var2name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_var3name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_var3name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_var3name as text
%        str2double(get(hObject,'String')) returns contents of edit_var3name as a double


% --- Executes during object creation, after setting all properties.
function edit_var3name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_var3name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_var4name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_var4name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_var4name as text
%        str2double(get(hObject,'String')) returns contents of edit_var4name as a double


% --- Executes during object creation, after setting all properties.
function edit_var4name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_var4name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_mainchoice.
function uipanel_mainchoice_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_mainchoice
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'togglebutton_math'
        handles.ImageIndex = 1;
        set(handles.radiobutton_mathoverlay,'Value',1);
        handles.overlayIndex = 1;
    case 'togglebutton_template'
        handles.ImageIndex = 2;
        set(handles.radiobutton_templateoverlay,'Value',1);
        handles.overlayIndex = 2;
    case 'togglebutton_ch1'
        handles.ImageIndex = 3;
        set(handles.radiobutton_ch1overlay,'Value',1);
        handles.overlayIndex = 3;
    case 'togglebutton_ch2'
        handles.ImageIndex = 4;
        set(handles.radiobutton_ch2overlay,'Value',1);
        handles.overlayIndex = 4;
    case 'togglebutton_ch3'
        handles.ImageIndex = 5;
        set(handles.radiobutton_ch3overlay,'Value',1);
        handles.overlayIndex = 5;
end
guidata(hObject, handles);
handles = guidata(hObject);

tp=str2num(get(handles.edit_currentFrame,'String'));
updatecurrentImage(handles,tp);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pushbutton_tonextframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_tonextframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_plotindividual.
function pushbutton_plotindividual_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plotindividual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
scell = str2num(get(handles.edit_cellNo,'String'));
allOutputs = get(handles.popupmenu_output,'String');
if ~strcmp(allOutputs,'<not present>')
    
    signal_name = ['/field' num2str(field)  '/' allOutputs{get(handles.popupmenu_output,'Value')}];
    outputNo = regexp(signal_name, ['(?<=.outputsignal)\d+'], 'match');
    if ~isempty(outputNo)
        timestamp_name = ['/field' num2str(field) '/timestamp' outputNo{1}];
    else
        timestamp_name = ['/field' num2str(field) '/timestamp'];
    end
    
    H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
    
    if exist(fullfile(handles.ndpathname,H5filename),'file')
        fileattrib(fullfile(handles.ndpathname,H5filename),'+w');
        fid = H5F.open(fullfile(handles.ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
        if H5L.exists(fid,signal_name,'H5P_DEFAULT') &&  H5L.exists(fid,timestamp_name,'H5P_DEFAULT')
            H5F.close(fid);
            legendList=[];
            signalinfo = h5info(fullfile(handles.ndpathname,H5filename), signal_name);
            startind = double([1 1 1]);
            countind = [signalinfo.Dataspace.Size(1) signalinfo.Dataspace.Size(2) 4];
            signal = permute(h5read(fullfile(handles.ndpathname,H5filename),signal_name,startind, countind),[2 1 3]);
            timestamp = double(h5read(fullfile(handles.ndpathname,H5filename),timestamp_name));
            
            figure(1);
            legendListCount = 1;
            
            if get(handles.checkbox_variable1,'Value')
                PosTime = find(signal(:,scell,1));
                if ~isempty(PosTime)
                    myy = signal(PosTime,scell,1)/mean(signal(PosTime,scell,1));
                    plot(timestamp(PosTime)/60,myy,'b');
                    if length(signalinfo.Attributes)~=0
                    legendList{legendListCount} = h5readatt(fullfile(handles.ndpathname,H5filename),signal_name,'signal1');
                    end
                    assignin('base','signal1',signal(PosTime,scell,1));
                    assignin('base','timestamp1',timestamp(PosTime));
                    legendListCount=legendListCount+1;
                else
                    display('Parameter 1 not available');
                end
            end
            
            if get(handles.checkbox_variable2,'Value')
                PosTime = find(signal(:,scell,2));
                if ~isempty(PosTime)
                    myy = signal(PosTime,scell,2)/mean(signal(PosTime,scell,2));
                    hold on;plot(timestamp(PosTime)/60,myy,'r');hold off;
                    if length(signalinfo.Attributes)~=0
                        legendList{legendListCount} = h5readatt(fullfile(handles.ndpathname,H5filename),signal_name,'signal2');
                    end
                    assignin('base','signal2',signal(PosTime,scell,2));
                    assignin('base','timestamp2',timestamp(PosTime));
                    legendListCount=legendListCount+1;
                else
                    display('Parameter 2 not available');
                end
                
            end
            if get(handles.checkbox_variable3,'Value')
                PosTime = find(signal(:,scell,3));
                if ~isempty(PosTime)
                    
                    myy = signal(PosTime,scell,3)/mean(signal(PosTime,scell,3));
                    hold on;plot(timestamp(PosTime)/60,myy,'g');hold off;
                    if length(signalinfo.Attributes)~=0
                        legendList{legendListCount} = h5readatt(fullfile(handles.ndpathname,H5filename),signal_name,'signal3');
                    end
                    assignin('base','signal3',signal(PosTime,scell,3));
                    assignin('base','timestamp3',timestamp(PosTime));
                    legendListCount=legendListCount+1;
                else
                    display('Parameter 3 not available');
                end
            end
            if get(handles.checkbox_variable4,'Value')
                PosTime = find(signal(:,scell,4));
                if ~isempty(PosTime)
                    myy = signal(PosTime,scell,4)/mean(signal(PosTime,scell,4));
                    hold on;plot(timestamp(PosTime)/60,myy,'k');hold off;
                    if length(signalinfo.Attributes)~=0
                        legendList{legendListCount} = h5readatt(fullfile(handles.ndpathname,H5filename),signal_name,'signal4');
                    end
                    assignin('base','signal4',signal(PosTime,scell,4));
                    assignin('base','timestamp4',timestamp(PosTime));
                    legendListCount=legendListCount+1;
                else
                    display('Parameter 4 not available');
                end
            end
            
            if get(handles.checkbox_variable1,'Value')
                
                myy = signal(:,scell,4)./signal(:,scell,3);
                myy(signal(:,scell,4) == 0 | signal(:,scell,3) == 0) = NaN;
                
                myy = myy/nanmean(myy);
                
                hold on;plot(timestamp/60,myy,'m');hold off;
                if length(signalinfo.Attributes)~=0
                    legendList{legendListCount} = [h5readatt(fullfile(handles.ndpathname,H5filename),signal_name,'signal4') '/' h5readatt(fullfile(handles.ndpathname,H5filename),signal_name,'signal3')];
                end
                assignin('base','signal5',myy);
                assignin('base','timestamp5',timestamp);
                legendListCount=legendListCount+1;
                
            end
            
            if length(signalinfo.Attributes)~=0
                title(h5readatt(fullfile(handles.ndpathname,H5filename),signal_name,'outputsignal_name'),'interpreter','none');
            end
            if legendListCount>1
                figure(1);
                legend(legendList);
            end
            xlabel('Time(hour)');
            
        else
            set(handles.edit_commu,'String',[signal_name 'or' timestamp_name ' does not exist']);
        end
    else
        set(handles.edit_commu,'String','Check to make sure that H5 file exists');
    end
else
    set(handles.edit_commu,'String','Data not available');
end

function text_ch1_Callback(hObject, eventdata, handles)
% hObject    handle to text_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_ch1 as text
%        str2double(get(hObject,'String')) returns contents of text_ch1 as a double
set(handles.togglebutton_ch1,'String',get(hObject,'String'));
mylist = get(handles.popupmenu_nomin,'String');
mylist{1} = get(hObject,'String');
set(handles.popupmenu_nomin,'String',mylist);

mylist = get(handles.popupmenu_denomin,'String');
mylist{1} = get(hObject,'String');
set(handles.popupmenu_denomin,'String',mylist);

mylist = get(handles.popupmenu_signal1,'String');
mylist{2} = get(hObject,'String');
set(handles.popupmenu_signal1,'String',mylist);
set(handles.popupmenu_signal2,'String',mylist);
set(handles.popupmenu_signal3,'String',mylist);
set(handles.popupmenu_signal4,'String',mylist);

function text_ch2_Callback(hObject, eventdata, handles)
% hObject    handle to text_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_ch2 as text
%        str2double(get(hObject,'String')) returns contents of text_ch2 as a double
set(handles.togglebutton_ch2,'String',get(hObject,'String'));
mylist = get(handles.popupmenu_nomin,'String');
mylist{2} = get(hObject,'String');
set(handles.popupmenu_nomin,'String',mylist);

mylist = get(handles.popupmenu_denomin,'String');
mylist{2} = get(hObject,'String');
set(handles.popupmenu_denomin,'String',mylist);

mylist = get(handles.popupmenu_signal1,'String');
mylist{3} = get(hObject,'String');
set(handles.popupmenu_signal1,'String',mylist);
set(handles.popupmenu_signal2,'String',mylist);
set(handles.popupmenu_signal3,'String',mylist);
set(handles.popupmenu_signal4,'String',mylist);

function text_ch3_Callback(hObject, eventdata, handles)
% hObject    handle to text_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_ch3 as text
%        str2double(get(hObject,'String')) returns contents of text_ch3 as a double
set(handles.togglebutton_ch3,'String',get(hObject,'String'));
mylist = get(handles.popupmenu_nomin,'String');
mylist{3} = get(hObject,'String');
set(handles.popupmenu_nomin,'String',mylist);

mylist = get(handles.popupmenu_denomin,'String');
mylist{3} = get(hObject,'String');
set(handles.popupmenu_denomin,'String',mylist);

mylist = get(handles.popupmenu_signal1,'String');
mylist{4} = get(hObject,'String');
set(handles.popupmenu_signal1,'String',mylist);
set(handles.popupmenu_signal2,'String',mylist);
set(handles.popupmenu_signal3,'String',mylist);
set(handles.popupmenu_signal4,'String',mylist);


% --- Executes when selected object is changed in uipanel_overlaychoice.
function uipanel_overlaychoice_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_overlaychoice
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton_mathoverlay'
        handles.overlayIndex = 1;
    case 'radiobutton_templateoverlay'
        
        handles.overlayIndex = 2;
    case 'radiobutton_ch1overlay'
        
        handles.overlayIndex = 3;
    case 'radiobutton_ch2overlay'
        
        handles.overlayIndex = 4;
    case 'radiobutton_ch3overlay'
        
        handles.overlayIndex = 5;
end
guidata(hObject, handles);
handles = guidata(hObject);

tp=str2num(get(handles.edit_currentFrame,'String'));
updatecurrentImage(handles,tp);
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

% --- Executes on button press in pushbutton_genmask.
function pushbutton_genmask_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_genmask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
first_tp = str2num(get(handles.edit_firstframe,'String'));
last_tp = str2num(get(handles.edit_lastframe,'String'));
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
nucCH= str2num(get(handles.edit_nucCH,'String'));
fileformat = get(handles.edit_signalformat,'String');
selected_cell = get(handles.listbox_cells,'Value');
CH1= handles.CH1;
CH2= handles.CH2;
CH3 = handles.CH3;
nucmask  = handles.nucmask;
cytomask = handles.cytomask;
cellmask = handles.cellmask;
stageLoc = get(handles.popupmenu_stagePosData,'Value');
ind_cellpath = handles.ind_cellpath;


switch get(handles.popupmenu_nomin,'Value')
    case 1
        nominCH=CH1;
    case 2
        nominCH=CH2;
    case 3
        nominCH=CH3;
end

switch get(handles.popupmenu_denomin,'Value')
    case 1
        denominCH=CH1;
    case 2
        denominCH=CH2;
    case 3
        denominCH=CH3;
    case 4
        denominCH=-1;
end

maskOUTfilename = ['maskOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(nucCH)];
load(maskOUTfilename,'selected_cells');

for tp=1:size(ind_cellpath,1)
    if ~isempty(find(selected_cell == selected_cells,1)) && ~isempty(find(nucmask(:,:,tp),1))
        NucMask  = nucmask(:,:,tp) ;
        
        if ~isempty(find(cellmask(:,:,tp),1))
            CellMask = cellmask(:,:,tp);
        else
            CellMask = bwmorph(NucMask,'dilate',str2num(get(handles.edit_cytosize,'String')));
            cellmask(:,:,tp) = CellMask;
        end
        
        if ~isempty(find(cytomask(:,:,tp),1))
            CytoMask = cytomask(:,:,tp);
        else
            CytoMask = bwmorph(CellMask,'erode',1)-NucMask;
            cytomask(:,:,tp) = CytoMask;
        end
    end
end

switch handles.filetype
    case 3
        MaskGenerating(3,[nucCH nominCH denominCH 1 last_tp],...
            [row col field plane],handles.channelnames, ...
            handles.ndfilename,selected_cell,stageLoc,...
            ind_cellpath,nucmask,cytomask,cellmask,handles.selected_cells);
end
handles.nucmask = nucmask;
handles.cellmask = cellmask;
handles.cytomask = cytomask;

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
            wavename3  = regexp(tline, '(?<="WaveName\d+", ")\w+(?=")', 'match');
            if ~isempty(wavename1) && ~isempty(wavename2)
                waveName{wind} = ['w' num2str(wind) wavename1{1} '-' wavename2{1}];
            else
                waveName{wind} = ['w' num2str(wind) wavename3{1}];
            end
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
[notp stagePos stageName waveName] = readndfile(handles.ndpathname,handles.ndfilename);

if notp==-1
    return;
end

set(handles.edit_firstframe,'String',num2str(1));
set(handles.edit_lastframe,'String',num2str(notp));
set(handles.edit_currentFrame,'String',num2str(1));
set(handles.popupmenu_stagePosData,'String',stagePos);
set(handles.popupmenu_stagePosData,'Value',1);
set(handles.popupmenu_stagePosBG,'String',stagePos);
set(handles.popupmenu_stagePosBG,'Value',1);
set(handles.edit_stageInfodata,'String',stageName{1});
set(handles.edit_stageInfobg,'String',stageName{1});
set(handles.edit_selectedSites,'String',['1:' num2str(length(stageName))]);

handles.stageName = stageName;
handles.channelnames = waveName;
set(handles.edit_totalCHs,'String',num2str(length(waveName)));
prefix = handles.ndfilename(1:(end-3));
handles.prefix = prefix;
fileformat = [prefix '_%s_s1_t%g.TIF'];
set(handles.edit_signalformat,'String',fileformat);

fileformat = [prefix '_%s_s1_t%g.TIF'];
set(handles.edit_fileformatBG,'String',fileformat);

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

handles.filetype = 3;
set(handles.radiobutton_customtiff,'Value',1);
handles.framestamp = 2;
set(handles.radiobutton_timestamp,'Value',1);
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_stagePosData.
function popupmenu_stagePosData_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_stagePosData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_stagePosData contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_stagePosData
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

switch handles.filetype
    case 3
        
        set(handles.edit_stageInfodata,'String',handles.stageName{get(hObject,'Value')});
        fileformat = [handles.prefix '_%s_s' num2str(get(hObject,'Value')) '_t%g.TIF'];
        set(handles.edit_signalformat,'String',fileformat);
        
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
        
    case 1
        row = str2num(get(handles.edit_row,'String'));
        col = str2num(get(handles.edit_col,'String'));
        field = str2num(get(handles.edit_field,'String'));
        plane = str2num(get(handles.edit_plane,'String'));
        set(handles.edit_signalformat,'String',['r' num2str(row,'%02.0f') 'c' num2str(col,'%02.0f') 'f' num2str(field,'%02.0f') 'p' num2str(plane,'%02.0f') 'rc%1.0f-ch1sk%ufk1fl1.tiff']);
end

guidata(hObject, handles);


% --- Executes on selection change in popupmenu_stagePosBG.
function popupmenu_stagePosBG_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_stagePosBG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_stagePosBG contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_stagePosBG
if handles.filetype==0
    set(handles.edit_commu,'String','No input images.');
    return;
end

switch handles.filetype
    case 3
        
        set(handles.edit_stageInfobg,'String',handles.stageName{get(hObject,'Value')});
        fileformat = [handles.prefix '_%s_s' num2str(get(hObject,'Value')) '_t%g.TIF'];
        set(handles.edit_fileformatBG,'String',fileformat);
        
    case 1
        row = str2num(get(handles.edit_row,'String'));
        col = str2num(get(handles.edit_col,'String'));
        field = str2num(get(handles.edit_field,'String'));
        plane = str2num(get(handles.edit_plane,'String'));
        set(handles.edit_fileformatBG,'String',['r' num2str(row,'%02.0f') 'c' num2str(col,'%02.0f') 'f' num2str(field,'%02.0f') 'p' num2str(plane,'%02.0f') 'rc%1.0f-ch1sk%ufk1fl1.tiff']);
end

guidata(hObject, handles);


%
% if handles.filetype==0
%     set(handles.edit_commu,'String','No input images.');
%     return;
% end
%
% set(handles.edit_stageInfobg,'String',handles.stageName{get(hObject,'Value')});
% fileformat = [handles.prefix '_%s_s' num2str(get(hObject,'Value')) '_t%g.TIF'];
% set(handles.edit_fileformatBG,'String',fileformat);
%
% filename = sprintf(fileformate,channelnames{channel},tp);
% h = fspecial('gaussian',[10 10],0.5);
% I = imfilter(imresize(im2double(imread('01152013-r2_w1Camera-YFP_s13_t1.TIF')),0.1),h,'replicate');
% [X,Y] = meshgrid(1:size(I,2),1:size(I,1));
% sfit = fit([X(:),Y(:)],I(:), 'poly33','Robust','LAR');
%
% guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_stagePosData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_stagePosData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stageInfodata_Callback(hObject, ~, handles)
% hObject    handle to edit_stageInfodata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stageInfodata as text
%        str2double(get(hObject,'String')) returns contents of edit_stageInfodata as a double


% --- Executes during object creation, after setting all properties.
function edit_stageInfodata_CreateFcn(hObject, ~, handles)
% hObject    handle to edit_stageInfodata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_squaresize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_squaresize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function uipanel_mainchoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_mainchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox_cellmarker.
function checkbox_cellmarker_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cellmarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_cellmarker

% --- Executes during object creation, after setting all properties.
function popupmenu_stagePosBG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_stagePosBG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_stageInfobg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stageInfobg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stageInfobg as text
%        str2double(get(hObject,'String')) returns contents of edit_stageInfobg as a double


% --- Executes during object creation, after setting all properties.
function edit_stageInfobg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stageInfobg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_fileformatBG_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fileformatBG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fileformatBG as text
%        str2double(get(hObject,'String')) returns contents of edit_fileformatBG as a double


% --- Executes during object creation, after setting all properties.
function edit_fileformatBG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fileformatBG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_useblank.
function checkbox_useblank_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_useblank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_useblank

if get(hObject,'Value')
    
    blankformat = get(handles.edit_fileformatBG,'String');
    row = str2num(get(handles.edit_row,'String'));
    col = str2num(get(handles.edit_col,'String'));
    field = str2num(get(handles.edit_field,'String'));
    plane = str2num(get(handles.edit_plane,'String'));
    channelnames = handles.channelnames;
    nucCH= str2num(get(handles.edit_nucCH,'String'));
    totalCHs = str2num(get(handles.edit_totalCHs,'String'));
    tp = str2num(get(handles.edit_currentFrame,'String'));
    blackIm = loadblank(handles,nucCH,tp);
    
    if ~isempty(blackIm)
        handles.useblank = 1;
        load fftexecutiontimes;
        handles.smooth_opt = detbestlength2(FFTrv,FFTiv,IFFTiv,size(blackIm),size(handles.gaussian),isreal(blackIm),isreal(handles.gaussian));
        set(handles.edit_commu,'String','Blank assigned');
        
    else
        
        handles.useblank = 0;
        set(handles.edit_commu,'String','Blank image does not exist');
    end
    
else
    handles.useblank = 0;
    set(handles.edit_commu,'String','Blank removed');
end

guidata(hObject, handles);

function edit_mark_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mark as text
%        str2double(get(hObject,'String')) returns contents of edit_mark as a double


% --- Executes during object creation, after setting all properties.
function edit_mark_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit64_Callback(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit64 as text
%        str2double(get(hObject,'String')) returns contents of edit64 as a double


% --- Executes during object creation, after setting all properties.
function edit64_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_autothres_math.
function pushbutton_autothres_math_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_autothres_math (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CH1= handles.CH1;
CH2= handles.CH2;
CH3 = handles.CH3;

tp = str2num(get(handles.edit_currentFrame,'String'));

switch get(handles.popupmenu_nomin,'Value')
    case 1
        nominCH=CH1;
    case 2
        nominCH=CH2;
    case 3
        nominCH=CH3;
end

switch get(handles.popupmenu_denomin,'Value')
    case 1
        denominCH=CH1;
    case 2
        denominCH=CH2;
    case 3
        denominCH=CH3;
    case 4
        denominCH=-1;
end
ratioIm = calculateFRET(handles,tp,nominCH,denominCH,handles.bg);
set(handles.edit_mathL,'String',num2str(min(ratioIm(:))));
set(handles.edit_mathH,'String',num2str(max(ratioIm(:))));

set(handles.togglebutton_math,'Value',1);
set(handles.radiobutton_mathoverlay,'Value',1);
handles.ImageIndex = 1;
handles.overlayIndex = 1;
guidata(hObject, handles);
handles = guidata(hObject);
updatecurrentImage(handles,tp);

function edit_templateL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_templateL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_templateL as text
%        str2double(get(hObject,'String')) returns contents of edit_templateL as a double


% --- Executes during object creation, after setting all properties.
function edit_templateL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_templateL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_templateH_Callback(hObject, eventdata, handles)
% hObject    handle to edit_templateH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_templateH as text
%        str2double(get(hObject,'String')) returns contents of edit_templateH as a double

% --- Executes during object creation, after setting all properties.
function edit_templateH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_templateH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_autothres_template.
function pushbutton_autothres_template_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_autothres_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tp = str2num(get(handles.edit_currentFrame,'String'));
channel = str2num(get(handles.edit_nucCH,'String'));
displayIM = loadsignal(handles,channel,tp);
set(handles.edit_templateL,'String',num2str(min(displayIM(:))));
set(handles.edit_templateH,'String',num2str(max(displayIM(:))));

set(handles.togglebutton_template,'Value',1);
set(handles.radiobutton_templateoverlay,'Value',1);
handles.ImageIndex = 2;
handles.overlayIndex = 2;
guidata(hObject, handles);
handles = guidata(hObject);
updatecurrentImage(handles,tp);

function edit_ch1L_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ch1L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ch1L as text
%        str2double(get(hObject,'String')) returns contents of edit_ch1L as a double


% --- Executes during object creation, after setting all properties.
function edit_ch1L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ch1L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ch1H_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ch1H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ch1H as text
%        str2double(get(hObject,'String')) returns contents of edit_ch1H as a double


% --- Executes during object creation, after setting all properties.
function edit_ch1H_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ch1H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_autothres_ch1.
function pushbutton_autothres_ch1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_autothres_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tp = str2num(get(handles.edit_currentFrame,'String'));
channel = handles.CH1;
displayIM = loadsignal(handles,channel,tp);
set(handles.edit_ch1L,'String',num2str(min(displayIM(:))));
set(handles.edit_ch1H,'String',num2str(max(displayIM(:))));
set(handles.togglebutton_ch1,'Value',1);
set(handles.radiobutton_ch1overlay,'Value',1);
handles.ImageIndex = 3;
handles.overlayIndex = 3;
guidata(hObject, handles);
handles = guidata(hObject);
updatecurrentImage(handles,tp);

function edit_ch2L_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ch2L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ch2L as text
%        str2double(get(hObject,'String')) returns contents of edit_ch2L as a double


% --- Executes during object creation, after setting all properties.
function edit_ch2L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ch2L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ch2H_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ch2H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ch2H as text
%        str2double(get(hObject,'String')) returns contents of edit_ch2H as a double


% --- Executes during object creation, after setting all properties.
function edit_ch2H_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ch2H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_autothres_ch2.
function pushbutton_autothres_ch2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_autothres_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tp = str2num(get(handles.edit_currentFrame,'String'));
channel = handles.CH2;
displayIM = loadsignal(handles,channel,tp);
set(handles.edit_ch2L,'String',num2str(min(displayIM(:))));
set(handles.edit_ch2H,'String',num2str(max(displayIM(:))));
set(handles.togglebutton_ch2,'Value',1);
set(handles.radiobutton_ch2overlay,'Value',1);
handles.ImageIndex = 4;
handles.overlayIndex = 4;
guidata(hObject, handles);
handles = guidata(hObject);
updatecurrentImage(handles,tp);

function edit_ch3L_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ch3L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ch3L as text
%        str2double(get(hObject,'String')) returns contents of edit_ch3L as a double


% --- Executes during object creation, after setting all properties.
function edit_ch3L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ch3L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ch3H_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ch3H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ch3H as text
%        str2double(get(hObject,'String')) returns contents of edit_ch3H as a double


% --- Executes during object creation, after setting all properties.
function edit_ch3H_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ch3H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_autothres_ch3.
function pushbutton_autothres_ch3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_autothres_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tp = str2num(get(handles.edit_currentFrame,'String'));
channel = handles.CH3;
displayIM = loadsignal(handles,channel,tp);
set(handles.edit_ch3L,'String',num2str(min(displayIM(:))));
set(handles.edit_ch3H,'String',num2str(max(displayIM(:))));
set(handles.togglebutton_ch3,'Value',1);
set(handles.radiobutton_ch3overlay,'Value',1);
handles.ImageIndex = 5;
handles.overlayIndex = 5;
guidata(hObject, handles);
handles = guidata(hObject);
updatecurrentImage(handles,tp);

function edit_BGsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_BGsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_BGsize as text
%        str2double(get(hObject,'String')) returns contents of edit_BGsize as a double


% --- Executes during object creation, after setting all properties.
function edit_BGsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_BGsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cytosize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cytosize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cytosize as text
%        str2double(get(hObject,'String')) returns contents of edit_cytosize as a double


% --- Executes on button press in checkbox_celllabeling.
function checkbox_celllabeling_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_celllabeling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_celllabeling


% --- Executes on button press in checkbox_cellnumber.
function checkbox_cellnumber_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cellnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_cellnumber



function edit_cellNo_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cellNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cellNo as text
%        str2double(get(hObject,'String')) returns contents of edit_cellNo as a double
set(handles.listbox_cells,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);



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
set(hObject,'String',pwd);
handles.ndpathname = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pushbutton_locatendfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_locatendfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in popupmenu_outputlist.
function popupmenu_outputlist_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_outputlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_outputlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_outputlist


% --- Executes during object creation, after setting all properties.
function popupmenu_outputlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_outputlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_removeoutputsignal.
function pushbutton_removeoutputsignal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removeoutputsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
allOutputs = get(handles.popupmenu_output,'String');
if ~strcmp(allOutputs,'<not present>')
    
    signal_name = ['/field' num2str(field)  '/' allOutputs{get(handles.popupmenu_output,'Value')}];
    outputNo = regexp(signal_name, ['(?<=.outputsignal)\d+'], 'match');
    if ~isempty(outputNo)
        timestamp_name = ['/field' num2str(field) '/timestamp' outputNo{1}];
    else
        timestamp_name = ['/field' num2str(field) '/timestamp'];
    end
    
    H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
    fileattrib(fullfile(handles.ndpathname,H5filename),'+w');
    fid = H5F.open(fullfile(handles.ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    
    if H5L.exists(fid,signal_name,'H5P_DEFAULT')
        H5L.delete(fid,signal_name,'H5P_DEFAULT');
        display(['Deleting output dataset ' H5filename ':' signal_name]);
    end
    if H5L.exists(fid,timestamp_name,'H5P_DEFAULT')
        H5L.delete(fid,timestamp_name,'H5P_DEFAULT');
        display(['Deleting output dataset ' H5filename ':' timestamp_name]);
    end
    
    H5F.close(fid);
    updateOutputList(handles);
end

% --- Executes on selection change in popupmenu_output.
function popupmenu_output_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_output contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_output
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
allOutputs = get(hObject,'String');
if ~strcmp(allOutputs,'<not present>')
    currentoutput_name = ['/field' num2str(field)  '/' allOutputs{get(hObject,'Value')}];
    H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
    
    myinfo=h5info(fullfile(handles.ndpathname,H5filename),currentoutput_name);
    %tmp = regexp(info.Datasets(i).Name, ['(?<=' basename ')\d+'], 'match');
    attLength = length(myinfo.Attributes);
    if attLength>0
        for i=1:attLength
            if strcmp(myinfo.Attributes(i).Name,'outputsignal_name')
                current_outname = myinfo.Attributes(i).Value;
                set(handles.edit_outputname,'String',current_outname);
            end
        end
    else
        set(handles.edit_outputname,'String','no name');
    end
end
% --- Executes during object creation, after setting all properties.
function popupmenu_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_storeall.
function pushbutton_storeall_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_storeall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CH1= handles.CH1;
CH2= handles.CH2;
CH3 = handles.CH3;
switch get(handles.popupmenu_nomin,'Value')
    case 1
        nominCH=CH1;
    case 2
        nominCH=CH2;
    case 3
        nominCH=CH3;
end
handles.nominCH= nominCH;

switch get(handles.popupmenu_denomin,'Value')
    case 1
        denominCH=CH1;
    case 2
        denominCH=CH2;
    case 3
        denominCH=CH3;
    case 4
        denominCH=-1;
end
handles.denominCH= denominCH;
handles.nucCH = str2num(get(handles.edit_nucCH,'String'));
handles.cellCH = str2num(get(handles.edit_cellCH,'String'));
handles.cellsize = str2num(get(handles.edit_cellsize,'String'));
handles.totalCHs = str2num(get(handles.edit_totalCHs,'String'));
handles.signalformat = get(handles.edit_signalformat,'String');
handles.blankformat = get(handles.edit_fileformatBG,'String');
handles.first_tp = str2num(get(handles.edit_firstframe,'String'));
handles.last_tp = str2num(get(handles.edit_lastframe,'String'));
handles.minstamp = get(handles.edit_minstamp,'String');
handles.secstamp = get(handles.edit_secstamp,'String');
handles.cytosize = get(handles.edit_cytosize,'String');

handles.bgnominType = get(handles.popupmenu_bgnomin,'Value');
handles.bgdenominType = get(handles.popupmenu_bgdenomin,'Value');
handles.illumlogic = get(handles.checkbox_illumlogic,'Value');

handles.Var1LOG = get(handles.checkbox_variable1,'Value');
handles.Var2LOG = get(handles.checkbox_variable2,'Value');
handles.Var3LOG = get(handles.checkbox_variable3,'Value');
handles.Var4LOG = get(handles.checkbox_variable4,'Value');

handles.region1LOG = get(handles.popupmenu_regionVar1,'Value');
handles.region2LOG = get(handles.popupmenu_regionVar2,'Value');
handles.region3LOG = get(handles.popupmenu_regionVar3,'Value');
handles.region4LOG = get(handles.popupmenu_regionVar4,'Value');

handles.signal1LOG = get(handles.popupmenu_signal1,'Value');
handles.signal2LOG = get(handles.popupmenu_signal2,'Value');
handles.signal3LOG = get(handles.popupmenu_signal3,'Value');
handles.signal4LOG = get(handles.popupmenu_signal4,'Value');

handles.regions1 = get(handles.popupmenu_regionVar1,'String');
handles.regions2 = get(handles.popupmenu_regionVar2,'String');
handles.regions3 = get(handles.popupmenu_regionVar3,'String');
handles.regions4 = get(handles.popupmenu_regionVar4,'String');

handles.signals1 = get(handles.popupmenu_signal1,'String');
handles.signals2 = get(handles.popupmenu_signal2,'String');
handles.signals3 = get(handles.popupmenu_signal3,'String');
handles.signals4 = get(handles.popupmenu_signal4,'String');

handles.outputsignaldatasetname = get(handles.edit_outputname,'String');

handles.filterParam1 = str2num(get(handles.edit_param1,'String'));
handles.filterParam2 = str2num(get(handles.edit_param2,'String'));
handles.bgsize = round(str2num(get(handles.edit_BGsize,'String'))/2);

handles.signalShiftN = str2num(get(handles.edit_bg_nomin_custom,'String'));
handles.signalShiftD = str2num(get(handles.edit_bg_denomin_custom,'String'));

guidata(hObject, handles);
handles = guidata(hObject);
BLANKsite = get(handles.popupmenu_stagePosBG,'Value');
[~, ~, stageName, ~] = readndfile(handles.ndpathname,handles.ndfilename);
sites = str2num(get(handles.edit_selectedSites,'String'));
if matlabpool('size') == 0
    matlabpool open;
end
parfor s=1:length(sites)
    site = sites(s);
    signalformat = [handles.prefix '_%s_s' num2str(site) '_t%g.TIF'];
    blankformat = [handles.prefix '_%s_s' num2str(BLANKsite) '_t%g.TIF'];
    L = regexp(stageName{site}, 'r(?<row>\d+)','names');
    if ~isempty(L)
        row = str2num(L.row);
    else
        row = site;
    end
    L = regexp(stageName{site}, 'c(?<col>\d+)','names');
    if ~isempty(L)
        col = str2num(L.col);
    else
        col = 1;
    end
    L = regexp(stageName{site}, 'f(?<field>\d+)','names');
    if ~isempty(L)
        field = str2num(L.field);
    else
        field = 1;
    end
    
    collectdata_individual(handles,row,col,field,signalformat,blankformat);
end
set(handles.edit_commu,'String','finished collecting all data');
matlabpool close;
updateOutputList(handles);

function collectdata_individual(handles,row,col,field,signalformat,blankformat)
nucCH = handles.nucCH;
CH1= handles.CH1;
CH2= handles.CH3;
CH3= handles.CH3;
outputsignaldatasetname = handles.outputsignaldatasetname;
nominCH= handles.nominCH;
denominCH= handles.denominCH;
first_tp=handles.first_tp ;
last_tp=handles.last_tp ;
totalCHs = handles.totalCHs;
cellsize = handles.cellsize;
Var1LOG=handles.Var1LOG ;
Var2LOG=handles.Var2LOG ;
Var3LOG=handles.Var3LOG ;
Var4LOG=handles.Var4LOG ;
region1LOG=handles.region1LOG ;
region2LOG=handles.region2LOG ;
region3LOG=handles.region3LOG;
region4LOG=handles.region4LOG;
signal1LOG=handles.signal1LOG ;
signal2LOG=handles.signal2LOG ;
signal3LOG=handles.signal3LOG ;
signal4LOG=handles.signal4LOG ;
filetype = handles.filetype;
channelnames = handles.channelnames;
ndpathname = handles.ndpathname;
outputsignalname=handles.outputsignalname;
nucCH = handles.nucCH; 
cellCH = handles.cellCH; 
H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
cellpath_name = ['/field' num2str(field) '/cellpath'];
sisterList_name = ['/field' num2str(field) '/sisterList'];
bg_name = ['/field' num2str(field) '/bg'];
maskdatasetname = ['/field' num2str(field) '/segmentsCH' num2str(nucCH)];
selectedcells_name = ['/field' num2str(field) '/selectedcells'];

if exist(fullfile(handles.ndpathname,H5filename),'file')
    fid = H5F.open(fullfile(handles.ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
        H5F.close(fid);
        cellpathinfo = h5info(fullfile(handles.ndpathname,H5filename), cellpath_name);
        cellpath_mat = h5read(fullfile(handles.ndpathname,H5filename),cellpath_name,[1 1 1], [cellpathinfo.Dataspace.Size(1) cellpathinfo.Dataspace.Size(2) cellpathinfo.Dataspace.Size(3)]);
        
        for tp=first_tp:last_tp
            cellpath{tp} = cellpath_mat(:,:,tp);
        end
    else
        H5F.close(fid);
        set(handles.edit_commu,'String',[H5filename ' ' cellpath_name 'does not exist.']);
        return;
    end
    
    fid = H5F.open(fullfile(handles.ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if H5L.exists(fid,maskdatasetname,'H5P_DEFAULT')
        H5F.close(fid);
        maskOption = 1;
    else
        nucFolder = [handles.ndpathname filesep 'nuclearMask'];
        cellFolder = [handles.ndpathname filesep 'cellMask'];
        if exist(nucFolder,'dir') && exist(cellFolder,'dir')
            
            maskOption = 2;
        else
            set(handles.edit_commu,'String','Mask folders do not exist.');
            return;
        end
    end
    
else
    set(handles.edit_commu,'String',[H5filename ' does not exist.']);
    return;
end

FractionExistingLimit = 0.4;
i=1;
serpen = [];
for c=1:size(cellpath_mat,1)
    myX = squeeze(cellpath_mat(c,1,:));
    myY = squeeze(cellpath_mat(c,2,:));
    pos_time = find(myX >0 & myY >0);
    
    if ~isempty(pos_time)
        
        end_p_idx = find(diff(pos_time)~=1,1,'first');
        start_p = min(pos_time);
        
        if isempty(end_p_idx)
            end_p = max(pos_time);
        else
            end_p = pos_time(end_p_idx);
        end
        serpen(i,:) = [end_p-start_p+1, c, start_p, end_p,];
        i=i+1;
    end
    clear myX myY pos_time;
end
serpen = sortrows(serpen,[-1 2 3 4]);
sSerpen_idx = find(serpen(:,1)> size(cellpath_mat,3)*FractionExistingLimit);
selected_cells = sort(serpen(sSerpen_idx,2));
fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,selectedcells_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' selectedcells_name]);
else
    H5L.delete(fid,selectedcells_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' selectedcells_name]);
    H5F.close(fid);
end

h5create(fullfile(ndpathname,H5filename), selectedcells_name, length(selected_cells), 'Datatype', 'uint32');
h5write(fullfile(ndpathname,H5filename), selectedcells_name, uint32(selected_cells(:)));


signal_name = ['/field' num2str(field) '/outputsignal1'];
timestamp_name = ['/field' num2str(field) '/timestamp1'];
basename='outputsignal';
info = h5info(fullfile(ndpathname,H5filename),['/field' num2str(field)]);
cindex = 1;
outputcounts=[];
if length(info.Datasets)>0 %#ok<*ISMT>
    for i=length(info.Datasets):-1:1
        outputsignalname = info.Datasets(i).Name;
        tmp = regexp(outputsignalname, ['(?<=' basename ')\d+'], 'match');
        if ~isempty(tmp)
            fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
            H5L.delete(fid,['/field' num2str(field) '/outputsignal' tmp{1}],'H5P_DEFAULT');
            H5F.close(fid);
            display(['Deleting r' num2str(row) 'c' num2str(col) 'f' num2str(field) '-' outputsignalname]);
        end
    end
end

% Initializing storage variables
mysignal = zeros(size(cellpath{last_tp},1),last_tp-first_tp+1,4);
h5create(fullfile(ndpathname,H5filename), signal_name, [size(mysignal,1), size(mysignal,2), 4], 'Datatype', 'double', 'ChunkSize', [size(mysignal,1), size(mysignal,2), 1], 'Deflate', 9);
h5write(fullfile(ndpathname,H5filename), signal_name, mysignal, [1 1 1], [size(mysignal,1) size(mysignal,2) 4]);
h5writeatt(fullfile(ndpathname,H5filename),signal_name,'outputsignal_name',outputsignaldatasetname);
if Var1LOG == 1
    region1 = handles.regions1{region1LOG};
    signal1 = handles.signals1{signal1LOG};
    h5writeatt(fullfile(ndpathname,H5filename),signal_name,'signal1',[region1 '-' signal1]);
end
if Var2LOG == 1
    region2 = handles.regions2{region2LOG};
    signal2 = handles.signals2{signal2LOG};
    h5writeatt(fullfile(ndpathname,H5filename),signal_name,'signal2',[region2 '-' signal2]);
end
if Var3LOG == 1
    region3 = handles.regions3{region3LOG};
    signal3 = handles.signals3{signal3LOG};
    h5writeatt(fullfile(ndpathname,H5filename),signal_name,'signal3',[region3 '-' signal3]);
end
if Var4LOG == 1
    region4 = handles.regions4{region4LOG};
    signal4 = handles.signals4{signal4LOG};
    h5writeatt(fullfile(ndpathname,H5filename),signal_name,'signal4',[region4 '-' signal4]);
end

timestamp = 1:(last_tp-first_tp+1);
fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,timestamp_name,'H5P_DEFAULT')
    H5F.close(fid);
else
    H5L.delete(fid,timestamp_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' timestamp_name]);
    H5F.close(fid);
end
h5create(fullfile(ndpathname,H5filename), timestamp_name, last_tp-first_tp+1, 'Datatype', 'double');

%[new_cellpath,new_sisterList] = removeSister(cellpath,sisterList,first_tp,last_tp,1:length(cellpath{last_tp}));

for tp=first_tp:last_tp
    % Determine image capture time based on the template channel
    if filetype == 3
        filename = sprintf(signalformat,channelnames{nucCH},first_tp);
        first_info = imfinfo(fullfile(ndpathname,filename));
        filename = sprintf(signalformat,channelnames{nucCH},tp);
        current_info = imfinfo(fullfile(ndpathname,filename));
    end
    
    switch handles.framestamp
        case 1
            timestep = tp-first_tp+1;
        case 2
            switch filetype
                case 3
                    [~, ~, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
                    hour = 24*D+H;
                    minute = MN;
                    second = round(S);
                otherwise
                    hour = floor(1/60*(tp-1)*( str2num(minstamp) + str2num(secstamp)/60 ));
                    minute = mod(floor((tp-1)*( str2num(minstamp) + str2num(secstamp)/60 )),60);
                    second = mod((tp-1)*( str2num(minstamp)*60 + str2num(secstamp)),60);
            end
            timestep = hour*60+ minute + second/60;
    end
    
    h5write(fullfile(ndpathname,H5filename), timestamp_name, timestep,tp-first_tp+1,1);
    % Determine signals
    
    skipframe = 0;
    switch maskOption
        case 1
            maskinfo = h5info(fullfile(ndpathname,H5filename), maskdatasetname);
        case 2
            if exist(fullfile(nucFolder, sprintf(signalformat,channelnames{nucCH},tp)),'file')
                nuc_im   = imread(fullfile(nucFolder, sprintf(signalformat,channelnames{nucCH},tp)));
                if exist(fullfile(cellFolder,sprintf(signalformat,channelnames{cellCH},tp)),'file')
                    cell_im  = imread(fullfile(cellFolder,sprintf(signalformat,channelnames{cellCH},tp)));
                elseif exist(fullfile(cellFolder,sprintf(signalformat,channelnames{nucCH},tp)),'file')
                    cell_im  = imread(fullfile(cellFolder,sprintf(signalformat,channelnames{nucCH},tp)));
                else
                    display([fullfile(cellFolder,sprintf(signalformat,channelnames{cellCH},tp)) ' does not exist.']);
                    skipframe = 1;
                end
            else
                display([fullfile(nucFolder,sprintf(signalformat,channelnames{cellCH},tp)) ' does not exist.']);
                skipframe = 1;
            end
    end
    if ~skipframe
        display(['Processing ' H5filename 'time point:' num2str(tp) ' of ' num2str(last_tp)]);
        CH1im = loadsignalV2(handles,CH1,tp,signalformat,blankformat);
        CH2im = loadsignalV2(handles,CH2,tp,signalformat,blankformat);
        CH3im = loadsignalV2(handles,CH3,tp,signalformat,blankformat);
        ratioIm = calculateFRETV2(handles,[],tp,nominCH,denominCH,signalformat,blankformat);
        
        imwidth = size(CH1im,2);
        imheight = size(CH1im,1);
        
        for s=sSerpen_idx'
            skipcell = 0;
            pos_time = (serpen(s,3):serpen(s,4));
            if ismember(tp,pos_time)
                cellNo = serpen(s,2);
                cX = cellpath_mat(cellNo,1,tp);
                cY = cellpath_mat(cellNo,2,tp);
                switch maskOption
                    case 1
                        startind = double([cellNo tp 1 1 1]);
                        countind = [1 1 3 maskinfo.Dataspace.Size(4) maskinfo.Dataspace.Size(5)];
                        allmasks = permute(h5read(fullfile(ndpathname,H5filename),maskdatasetname,startind, countind),[4 5 3 1 2]);
                        nucmask  = double(allmasks(:,:,1));
                        cellmask  = double(allmasks(:,:,2));
                        cytomask  = double(allmasks(:,:,3));
                        if isempty(find(nucmask,1))
                            skipcell = 1;
                        end
                        
                    case 2
                        cellmaskT = bwselect(cell_im,cX,cY,8);
                        nucmaskT  = bwselect(nuc_im,cX,cY,8);
                        if ~isempty(find(nucmaskT,1)) || ~isempty(find(cellmaskT,1))
                            S = regionprops(cellmaskT,{'BoundingBox'});
                            BoundingBox = round(S(1).BoundingBox);
                            
                            xL=max(BoundingBox(1),1);
                            xR=min(xL+BoundingBox(3),imwidth);
                            yL=max(BoundingBox(2),1);
                            yR=min(yL+BoundingBox(4),imheight);

                            cellmask = cellmaskT(yL:yR,xL:xR);
                            nucmask  = nucmaskT(yL:yR,xL:xR);
                            cytomask = cellmask-nucmask;
                        else
                            skipcell = 1;
                        end
                end
                if ~skipcell
                    %for Nuclei region
                    [nuc_Y, nuc_X] = find(nucmask~=0);
                    %for Cytosol region
                    [cyto_Y, cyto_X] = find(cytomask~=0);
                    %for Cell region
                    [cell_Y, cell_X] = find(cellmask~=0);

                    switch maskOption
                        case 1
                            if cellsize~=(size(nucmask,1)-1)/2
                                cellsize = (size(nucmask,1)-1)/2;
                            end
                            xL=max(cX-cellsize,1);
                            xR=min(cX+cellsize,imwidth);
                            yL=max(cY-cellsize,1);
                            yR=min(cY+cellsize,imheight);
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
                            cell_ch1 = zeros(2*cellsize+1,2*cellsize+1);
                            cell_ch2 = zeros(2*cellsize+1,2*cellsize+1);
                            cell_ch3 = zeros(2*cellsize+1,2*cellsize+1);
                            mini_ratioIm = zeros(2*cellsize+1,2*cellsize+1);
                            
                            cell_ch1(borderY,borderX) = CH1im(yL:yR,xL:xR);
                            cell_ch2(borderY,borderX) = CH2im(yL:yR,xL:xR);
                            cell_ch3(borderY,borderX) = CH3im(yL:yR,xL:xR);
                            mini_ratioIm(borderY,borderX) = ratioIm(yL:yR,xL:xR);
                        case 2
                            cell_ch1 = CH1im(yL:yR,xL:xR);
                            cell_ch2 = CH2im(yL:yR,xL:xR);
                            cell_ch3 = CH3im(yL:yR,xL:xR);
                            mini_ratioIm = ratioIm(yL:yR,xL:xR);
                    end
                    
                    if Var1LOG == 1
                        mysignal(cellNo,tp-first_tp+1,1) = signalOutputing(region1LOG,signal1LOG,mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y);
                    end
                    
                    if Var2LOG == 1
                        mysignal(cellNo,tp-first_tp+1,2) = signalOutputing(region2LOG,signal2LOG,mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y);
                    end
                    
                    if Var3LOG == 1
                        mysignal(cellNo,tp-first_tp+1,3) = signalOutputing(region3LOG,signal3LOG,mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y);
                    end
                    
                    if Var4LOG == 1
                        mysignal(cellNo,tp-first_tp+1,4) = signalOutputing(region4LOG,signal4LOG,mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y);
                    end
                    clear nuc_X nuc_Y cyto_X cyto_Y cell_X cell_Y cell_template mini_ratioIm cell_ch1 cell_ch2 cell_ch3
                    %else
                    %    display(['Skipping cell:' num2str(cellNo) H5filename 'time point:' num2str(tp) ]);
                end
            end
        end
        h5write(fullfile(ndpathname,H5filename), signal_name, mysignal(:,tp-first_tp+1,:), [1 tp-first_tp+1 1], [size(mysignal,1) 1 4]);
        
        %else
        %    display(['Skipping all cells ' H5filename 'time point:' num2str(tp) ' of ' num2str(last_tp)]);
        clear  ratioIm CH1im CH2im CH3im;
    end
end

display(['r' num2str(row) 'c' num2str(col) 'f' num2str(field) '-Finished']);



function edit_selectedSites_Callback(hObject, eventdata, handles)
% hObject    handle to edit_selectedSites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_selectedSites as text
%        str2double(get(hObject,'String')) returns contents of edit_selectedSites as a double


% --- Executes during object creation, after setting all properties.
function edit_selectedSites_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_selectedSites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mathname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mathname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mathname as text
%        str2double(get(hObject,'String')) returns contents of edit_mathname as a double

mylist = get(handles.popupmenu_signal1,'String');
mylist{1} = get(hObject,'String');
set(handles.popupmenu_signal1,'String',mylist);
set(handles.popupmenu_signal2,'String',mylist);
set(handles.popupmenu_signal3,'String',mylist);
set(handles.popupmenu_signal4,'String',mylist);

% --- Executes during object creation, after setting all properties.
function edit_mathname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mathname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cellCH_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cellCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cellCH as text
%        str2double(get(hObject,'String')) returns contents of edit_cellCH as a double


% --- Executes during object creation, after setting all properties.
function edit_cellCH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cellCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
