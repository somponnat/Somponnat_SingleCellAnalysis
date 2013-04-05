function varargout = FRETPlotting(varargin)
% Edit the above text to modify the response to help mytab

% Last Modified by GUIDE v2.5 23-Jan-2013 17:42:46

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
            
            if exist(signal_filename,'file') 
                signalim = im2double(imread(signal_filename));
                if exist(blank_filename,'file')
                    blankim = im2double(imread(blank_filename));
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
            if exist(signalformat,'file')
                outputim = im2double(imread(signalformat,'Index',totalCHs*(tp-1)+channel));
            else
                set(handles.edit_commu,'String',[signalformat ' does not exist.']);
            end
        case 3
            signal_filename = sprintf(signalformat,channelnames{channel},tp);
            blank_filename = sprintf(blankformat,channelnames{channel},tp);
            if exist(signal_filename,'file') 
                signalim = im2double(imread(signal_filename));
                if exist(blank_filename,'file')
                    blankim = im2double(imread(blank_filename));
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
            if exist(filename,'file')
                outputim = im2double(imread(filename));
            else
                set(handles.edit_commu,'String',[filename ' does not exist.']);
            end
        case 2
            if exist(signalformat,'file')
                outputim = im2double(imread(signalformat,'Index',totalCHs*(tp-1)+channel));
            else
                set(handles.edit_commu,'String',[signalformat ' does not exist.']);
            end
        case 3
            filename = sprintf(signalformat,channelnames{channel},tp);
            if exist(filename,'file');
                outputim = im2double(imread(filename));
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
        if exist(filename,'file')
            outputim = im2double(imread(filename));
        else
            set(handles.edit_commu,'String',[filename ' does not exist.']);
        end
    case 2
        if exist(signalformat,'file')
            outputim = im2double(imread(blankformat,'Index',totalCHs*(tp-1)+channel));
        else
            set(handles.edit_commu,'String',[signalformat ' does not exist.']);
        end
    case 3
        filename = sprintf(blankformat,channelnames{channel},tp);
        if exist(filename,'file');
            outputim = im2double(imread(filename));
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
        
        templateCH = trackinginfo(1);
        set(handles.edit_template,'String',num2str(templateCH));
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
        
        
        templateCH = trackinginfo(1);
        set(handles.edit_template,'String',num2str(templateCH));
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
            templateCH = 1;
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
            templateCH = trackinginfo(1);
            set(handles.edit_template,'String',num2str(templateCH));
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
set(handles.edit_celltrackfilename,'String','<cell track file>');
handles.celltrackpathname = [];
handles.fcn1 = makeConstrainToRectFcn('impoint',get(handles.axes1,'XLim'),get(handles.axes1,'YLim'));
handles.useblank = 0;


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
    
    template = loadsignal(handles,templateCH,tp);
    load fftexecutiontimes;
    h= 1/273*[1 4 7 4 1;4 16 26 26 4;7 26 41 26 7;4 16 26 16 4;1 4 7 4 1];
    handles.smooth_opt = detbestlength2(FFTrv,FFTiv,IFFTiv,size(template),size(h),isreal(template),isreal(h));
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
templateCH= handles.templateCH;
CH1= handles.CH1;
CH2= handles.CH2;
CH3 = handles.CH3;
cellsize = str2num(get(handles.edit_cellsize,'String')); %#ok<*NASGU>
tp=str2num(get(handles.edit_firstframe,'String')); %#ok<*ST2NM>
set(handles.edit_currentFrame,'String',num2str(tp));
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

filetype = handles.filetype;
channelnames = handles.channelnames;
template = loadsignal(handles,templateCH,tp);

startImageIndex = 1;
handles.ImageIndex = startImageIndex;
handles.overlayIndex = startImageIndex;
set(handles.togglebutton_math,'Value',1);
set(handles.radiobutton_mathoverlay,'Value',1);

fftw('planner', 'hybrid');
load MyColormaps;
handles.mycmap1 = mycmap1;
handles.mycmap2 = mycmap2;
handles.mycmap3 = mycmap3;
handles.mycmap4 = mycmap4;

switch startImageIndex
    case 3
        displayIM = loadsignal(handles,CH1,tp);
        if ~isempty(displayIM)
            imshow(displayIM,[],'Parent',handles.axes1);colormap gray;drawnow;
        end
    case 4
        displayIM = loadsignal(handles,CH2,tp);
        if ~isempty(displayIM)
            imshow(displayIM,[],'Parent',handles.axes1);colormap gray;drawnow;
        end
    case 5
        displayIM = loadsignal(handles,CH3,tp);
        if ~isempty(displayIM)
            imshow(displayIM,[],'Parent',handles.axes1);colormap gray;drawnow;
        end
    case 2
        displayIM = template;
        if ~isempty(displayIM)
            imshow(displayIM,[],'Parent',handles.axes1);colormap gray;drawnow;
        end
    case 1
        displayIM = calculateFRET(handles,tp,nominCH,denominCH,[]);
        if ~isempty(displayIM)
            displayratioIm(displayIM,handles);
        end
end

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
set(handles.edit_celltrackfilename,'String','<cell track file>');
handles.celltrackpathname = [];
handles.fcn1 = makeConstrainToRectFcn('impoint',get(handles.axes1,'XLim'),get(handles.axes1,'YLim'));
handles.useblank = 0;
set(handles.checkbox_useblank,'Value',0);
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
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
templateCH= handles.templateCH;
maskOUTfilename = ['maskOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(templateCH)];

nucName = ['nucmask_cell' num2str(selected_cell)];
cytoName = ['cytomask_cell' num2str(selected_cell)];
cellName = ['cellmask_cell' num2str(selected_cell)];

load(maskOUTfilename,nucName);
load(maskOUTfilename,cytoName);
load(maskOUTfilename,cellName);
load(maskOUTfilename,'selected_cells');

set(handles.edit_collectedcells,'String',num2str(selected_cells'));
set(handles.edit_commu,'String',[maskOUTfilename ': ' num2str(length(selected_cells)) ' cells selected']);

eval(['nucmask=' nucName ';']);
eval(['cytomask=' cytoName ';']);
eval(['cellmask=' cellName ';']);

clear(nucName);
clear(cytoName);
clear(cellName);

cellpath = handles.cellpath;
sisterList = handles.sisterList;
cellsize = str2num(get(handles.edit_cellsize,'String'));

tp = str2num(get(handles.edit_currentFrame,'String'));
firsttp = str2num(get(handles.edit_firstframe,'String'));
lasttp = str2num(get(handles.edit_lastframe,'String'));
currentframe = loadsignal(handles,templateCH,tp);
imwidth = size(currentframe,2);
imheight = size(currentframe,1);


ind_cellpath = pos_path(cellpath,sisterList,selected_cell,firsttp,lasttp,imheight,imwidth);
handles.ind_cellpath = ind_cellpath;

if tp > size(ind_cellpath,1)
    set(handles.edit_commu,'String',['Cell#' num2str(selected_cell) ' is not present']);
else

    if cellsize~=(size(nucmask(:,:,tp),1)-1)/2
        cellsize = (size(nucmask(:,:,tp),1)-1)/2;
    end
    
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
        bwCytoEdge = bwperim(CellMask);
        bwCellEdge = bwperim(CytoMask);
        
        bwFinal = bwNucEdge| bwCellEdge | bwCytoEdge;
        imAdj = imAdj-bwFinal;
        imOut=cat(3,max(imAdj,(bwNucEdge | cp)),max(imAdj,bwCellEdge),max(imAdj,bwCytoEdge));
    else
        set(handles.edit_commu,'String',['Cell#' num2str(selected_cell) ' has no masks']);
    end
    imshow(imOut,[],'Parent',handles.axes2);
    
    set(handles.edit_coord,'String',(sprintf('(%1.0f,%1.0f)',ind_cellpath(tp,1),ind_cellpath(tp,2))));
    set(handles.edit_cellNo,'String',num2str(selected_cell));
    set(handles.edit_sister,'String',num2str(sisterList{tp}(selected_cell,:)));

end
handles.nucmask = nucmask;
handles.cellmask = cellmask;
handles.cytomask = cytomask;

guidata(hObject, handles);

function  out_cellpath = pos_path(cellpath,sisterList,cellNo,firsttp,lasttp,imheight,imwidth)
% check if this path is good
%#1 if acquiring and being inferior sister, combine trace of self with
%prior-sister
% has sister?
new_cellpath = zeros(lasttp,2);

if  sisterList{end}(cellNo,1) ~= -1
    sis1No = sisterList{lasttp}(cellNo,1);
    
    if cellpath{firsttp}(cellNo,1) == -1
        mainInd = sis1No;
        subInd =  cellNo;
        subData = zeros(1,lasttp);
        for t=firsttp:lasttp
            subData(t) = sisterList{t}(subInd,1);
        end
        t_break = find(subData(firsttp:lasttp)>0,1,'first');
        
        for t=firsttp:lasttp
            if t<t_break
                new_cellpath(t,:) = cellpath{t}(mainInd,:);
            else
                new_cellpath(t,:) = cellpath{t}(subInd,:);
            end
        end
    else
        for t=firsttp:lasttp
            new_cellpath(t,:) = cellpath{t}(cellNo,:);
        end
    end
else
    for t=firsttp:lasttp
        new_cellpath(t,:) = cellpath{t}(cellNo,:);
    end
end

deathInd = find(new_cellpath(:,1)==-2,1,'first');

if isempty(deathInd)
    out_cellpath = new_cellpath;
else
    out_cellpath = new_cellpath(1:(deathInd-1),:);
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
    templateCH= handles.templateCH;
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
        filename = sprintf(fileformat,channelnames{templateCH},str2num(get(handles.edit_firstframe,'String')));
        first_info = imfinfo(filename);
        filename = sprintf(fileformat,channelnames{templateCH},tp);
        current_info = imfinfo(filename);
        [~, ~, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
        hour = 24*D+round(H);
        minute = round(MN);
        second = round(S);
        set(handles.edit_commu,'String',[num2str(hour,'%02.0f') ':' num2str(minute,'%02.0f') ':' num2str(second,'%02.0f')]);
    end
    
    filetype = handles.filetype;
    channelnames = handles.channelnames;
    template = loadsignal(handles,templateCH,tp);
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
                BG_N(b) = median(selectedN(:));
            end
        end
    case 3
        xL=round(size(nominIM,2)/2)-2*bgsize;
        xR=round(size(nominIM,2)/2)+2*bgsize;
        yL=round(size(nominIM,1)/2)-2*bgsize;
        yR=round(size(nominIM,1)/2)+2*bgsize;
        BG_N = median(nominBLK(:));
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
                BG_D(b) = median(selectedD(:));
            end
        end
    case 3
        xL=round(size(nominIM,2)/2)-2*bgsize;
        xR=round(size(nominIM,2)/2)+2*bgsize;
        yL=round(size(nominIM,1)/2)-2*bgsize;
        yR=round(size(nominIM,1)/2)+2*bgsize;
        BG_D = median(denomBLK(:));
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

load(fullfile(handles.celltrackpathname,get(handles.edit_celltrackfilename,'String')));

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

function edit_template_Callback(hObject, eventdata, handles)
% hObject    handle to edit_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_template as text
%        str2double(get(hObject,'String')) returns contents of edit_template as a double
handles.templateCH= str2num(get(hObject,'String'));
guidata(hObject, handles);  

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
handles.templateCH= str2num(get(hObject,'String'));
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
templateCH= handles.templateCH;
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
    filename = sprintf(fileformat,handles.channelnames{templateCH},first_tp);
    first_info = imfinfo(filename);
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
template = loadsignal(handles,templateCH,tp);
if handles.filetype == 3
    filename = sprintf(fileformat,handles.channelnames{templateCH},tp);
    current_info = imfinfo(filename);
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
templateCH= handles.templateCH;
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

aviobj = VideoWriter(['myMov_r' num2str(row) 'c' num2str(col) 'f' num2str(field) 'ch' num2str(handles.ImageIndex) '.avi']);
aviobj.FrameRate = str2num(get(handles.edit_framerate,'String'));
aviobj.Quality = 80;
open(aviobj);
f1=figure();
axes(gca);
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

if handles.filetype == 3
    filename = sprintf(fileformat,handles.channelnames{templateCH},first_tp);
    first_info = imfinfo(filename);
end


for tp=first_tp:last_tp
    
    
    if handles.ImageIndex == handles.overlayIndex
        
        template = loadsignal(handles,templateCH,tp);
        if handles.filetype == 3
            filename = sprintf(fileformat,handles.channelnames{templateCH},tp);
            current_info = imfinfo(filename);
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
        
        template = loadsignal(handles,templateCH,tp);
        if handles.filetype == 3
            filename = sprintf(fileformat,handles.channelnames{templateCH},tp);
            current_info = imfinfo(filename);
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
                I2 = mat2gray(loadsignal(handles,CH1,tp));
                
            case 4
                I2 = mat2gray(loadsignal(handles,CH2,tp));
                
            case 5
                I2 = mat2gray(loadsignal(handles,CH3,tp));
                
            case 2
                I2 = mat2gray(template);
                
            case 1
                I2 = mat2gray(calculateFRET(handles,tp,nominCH,denominCH,bg));
        end
        
        combinedI(:,:,1) = 0.8*I1;
        combinedI(:,:,2) = 0.8*I1+0.6*I2;
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

function edit_celltrackfilename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_celltrackfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_celltrackfilename as text
%        str2double(get(hObject,'String')) returns contents of edit_celltrackfilename as a double


% --- Executes during object creation, after setting all properties.
function edit_celltrackfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_celltrackfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_celltracklocator.
function pushbutton_celltracklocator_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_celltracklocator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname,FilterIndex] = uigetfile(...
{'*.mat','cell track mat file(*.mat)';
   '*.*',  'All Files (*.*)'}, 'Pick a file');
if FilterIndex~=0
    set(handles.edit_celltrackfilename,'String',filename);
    handles.celltrackpathname = pathname;
else
    set(handles.edit_celltrackfilename,'String','N/A');
    handles.celltrackpathname = [];
end
guidata(hObject, handles);  


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
templateCH= handles.templateCH;
CH1= handles.CH1;
CH2= handles.CH2;
CH3 = handles.CH3;

cellsize = str2num(get(handles.edit_cellsize,'String'));
totalCHs = str2num(get(handles.edit_totalCHs,'String'));
fileformat = get(handles.edit_signalformat,'String');
filetype = handles.filetype;
channelnames = handles.channelnames;
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

% Initializing storage variables
mysignal1 = cell(size(cellpath{last_tp},1),1);
mysignal2 = cell(size(cellpath{last_tp},1),1);
mysignal3 = cell(size(cellpath{last_tp},1),1);
mysignal4 = cell(size(cellpath{last_tp},1),1);

maskOUTfilename = ['maskOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(templateCH)];
load(maskOUTfilename,'selected_cells');
selected_cells = sort(selected_cells');
groupNo = 75;
for i=1:ceil(length(selected_cells)/groupNo)
    first_cell = groupNo*(i-1) + 1;
    if i < ceil(length(selected_cells)/groupNo)
        last_cell = first_cell+groupNo-1;
    else
        last_cell = length(selected_cells);
    end
    cell_list{i} = selected_cells(first_cell:last_cell);
end

for c_set = 1:length(cell_list)
    nucmask = cell(groupNo,1);
    cytomask = cell(groupNo,1);
    cellmask = cell(groupNo,1);
    current_list = cell_list{c_set};
    
    for loop_c=current_list
        nucName = ['nucmask_cell' num2str(loop_c)];
        cytoName = ['cytomask_cell' num2str(loop_c)];
        cellName = ['cellmask_cell' num2str(loop_c)];
        
        load(maskOUTfilename,nucName);
        load(maskOUTfilename,cytoName);
        load(maskOUTfilename,cellName);
        
        current_ind = mod(loop_c,groupNo);
        if current_ind== 0 
            current_ind = groupNo;
        end
        
        eval(['nucmask{' num2str(current_ind) '}=' nucName ';']);
        eval(['cytomask{' num2str(current_ind) '}=' cytoName ';']);
        eval(['cellmask{' num2str(current_ind) '}=' cellName ';']);
        
        clear(nucName);
        clear(cytoName);
        clear(cellName);
    end
  
    for tp=first_tp:last_tp
        % Determine image capture time based on the template channel
        if handles.filetype == 3
            filename = sprintf(fileformat,channelnames{templateCH},first_tp);
            first_info = imfinfo(filename);
            filename = sprintf(fileformat,channelnames{templateCH},tp);
            current_info = imfinfo(filename);
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
                        hour = floor(1/60*(tp-1)*( str2num(get(handles.edit_minstamp,'String')) + str2num(get(handles.edit_secstamp,'String'))/60 ));
                        minute = mod(floor((tp-1)*( str2num(get(handles.edit_minstamp,'String')) + str2num(get(handles.edit_secstamp,'String'))/60 )),60);
                        second = mod((tp-1)*( str2num(get(handles.edit_minstamp,'String'))*60 + str2num(get(handles.edit_secstamp,'String'))),60);
                end
                timestep = hour*60+ minute + second/60;
        end
        clc;display(['Set:' num2str(c_set) ' of ' num2str(length(cell_list)) ' - Processing time point:' num2str(tp) ' of ' num2str(last_tp)]);
        guidata(hObject, handles);
        handles = guidata(hObject);
        
        
        
        % Determine signals
        CH1im = loadsignal(handles,CH1,tp);
        CH2im = loadsignal(handles,CH2,tp);
        CH3im = loadsignal(handles,CH3,tp);
        ratioIm = calculateFRET(handles,tp,nominCH,denominCH,bg);
        
        imwidth = size(CH1im,2);
        imheight = size(CH1im,1);
        
        for cellNo=current_list
            
            ind_cellpath = pos_path(cellpath,sisterList,cellNo,first_tp,last_tp,imheight,imwidth);
            
            if tp <= size(ind_cellpath,1)
                
                if ~isempty(nucmask{current_ind}) && ~isempty(find(nucmask{current_ind}(:,:,tp),1))
                    
                    if cellsize~=(size(nucmask{current_ind}(:,:,first_tp),1)-1)/2
                        cellsize = (size(nucmask{current_ind}(:,:,first_tp),1)-1)/2;
                    end
                    
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
                    cell_ch1 = zeros(2*cellsize+1,2*cellsize+1);
                    cell_ch2 = zeros(2*cellsize+1,2*cellsize+1);
                    cell_ch3 = zeros(2*cellsize+1,2*cellsize+1);
                    mini_ratioIm = zeros(2*cellsize+1,2*cellsize+1);
                    
                    cell_ch1(borderY,borderX) = CH1im(yL:yR,xL:xR);
                    cell_ch2(borderY,borderX) = CH2im(yL:yR,xL:xR);
                    cell_ch3(borderY,borderX) = CH3im(yL:yR,xL:xR);
                    mini_ratioIm(borderY,borderX) = ratioIm(yL:yR,xL:xR);
                    
                    NucMask  = nucmask{current_ind}(:,:,tp) ;
                    
                    if ~isempty(find(cellmask{current_ind}(:,:,tp),1))
                        CellMask = cellmask{current_ind}(:,:,tp);
                    else
                        CellMask = bwmorph(NucMask,'dilate',str2num(get(handles.edit_cytosize,'String')));
                        cellmask{current_ind}(:,:,tp) = CellMask;
                    end
                    
                    if ~isempty(find(cytomask{current_ind}(:,:,tp),1))
                        CytoMask = cytomask{current_ind}(:,:,tp);
                    else
                        CytoMask = CellMask-NucMask;
                        cytomask{current_ind}(:,:,tp) = CytoMask;
                    end
                    
                    %for Nuclei region
                    [nuc_Y, nuc_X] = find(NucMask);
                    %for Cytosol region
                    [cyto_Y, cyto_X] = find(CytoMask);
                    %for Cell region
                    [cell_Y, cell_X] = find(CellMask);
                    
                    clear NucMask CytoMask CellMask;
                    
                    if get(handles.checkbox_variable1,'Value') == 1
                        mysignal1{cellNo} = [mysignal1{cellNo};timestep signalOutputing(get(handles.popupmenu_regionVar1,'Value'),get(handles.popupmenu_signal1,'Value'),mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y)];
                    end
                    
                    if get(handles.checkbox_variable2,'Value') == 1
                        mysignal2{cellNo} = [mysignal2{cellNo};timestep signalOutputing(get(handles.popupmenu_regionVar2,'Value'),get(handles.popupmenu_signal2,'Value'),mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y)];
                    end
                    
                    if get(handles.checkbox_variable3,'Value') == 1
                        mysignal3{cellNo} = [mysignal3{cellNo};timestep signalOutputing(get(handles.popupmenu_regionVar3,'Value'),get(handles.popupmenu_signal3,'Value'),mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y)];
                    end
                    
                    if get(handles.checkbox_variable4,'Value') == 1
                        mysignal4{cellNo} = [mysignal4{cellNo};timestep signalOutputing(get(handles.popupmenu_regionVar4,'Value'),get(handles.popupmenu_signal4,'Value'),mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y)];
                    end
                    clear nuc_X nuc_Y cyto_X cyto_Y cell_X cell_Y cell_template mini_ratioIm cell_ch1 cell_ch2 cell_ch3
                else
                    if get(handles.checkbox_variable1,'Value') == 1
                        mysignal1{cellNo} = [mysignal1{cellNo};timestep 0];
                    end
                    
                    if get(handles.checkbox_variable2,'Value') == 1
                        mysignal2{cellNo} = [mysignal2{cellNo};timestep 0];
                    end
                    
                    if get(handles.checkbox_variable3,'Value') == 1
                        mysignal3{cellNo} = [mysignal3{cellNo};timestep 0];
                    end
                    
                    if get(handles.checkbox_variable4,'Value') == 1
                        mysignal4{cellNo} = [mysignal4{cellNo};timestep 0];
                    end
                end
            end
        end
        
        clear  ratioIm CH1im CH2im CH3im
    end
end

savefile = ['AllFRETsignals_r' num2str(row) 'c' num2str(col) 'f' num2str(field)];
save(savefile,'selected_cells','-v7.3');

if get(handles.checkbox_variable1,'Value') == 1
    save(savefile,'mysignal1','-v7.3','-append');
end

if get(handles.checkbox_variable2,'Value') == 1
    save(savefile,'mysignal2','-v7.3','-append');
end

if get(handles.checkbox_variable3,'Value') == 1
    save(savefile,'mysignal3','-v7.3','-append');
end

if get(handles.checkbox_variable4,'Value') == 1
    save(savefile,'mysignal4','-v7.3','-append');
end
set(handles.edit_commu,'String',['All results saved in ' savefile '.mat']);
guidata(hObject, handles);

function outputsignal = signalOutputing(regiontype,signaltype,mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y)
switch regiontype
    case 1
        switch signaltype
            case 1
                Average_Nuc = mean(mean(mini_ratioIm(nuc_Y,nuc_X)));
                Average_Cyto = mean(mean(mini_ratioIm(cyto_Y,cyto_X)));

            case 2
                Average_Nuc = mean(mean(cell_ch1(nuc_Y,nuc_X)));
                Average_Cyto = mean(mean(cell_ch1(cyto_Y,cyto_X)));
                
            case 3
                Average_Nuc = mean(mean(cell_ch2(nuc_Y,nuc_X)));
                Average_Cyto = mean(mean(cell_ch2(cyto_Y,cyto_X)));
            case 4
                Average_Nuc = mean(mean(cell_ch3(nuc_Y,nuc_X)));
                Average_Cyto = mean(mean(cell_ch3(cyto_Y,cyto_X)));
        end
        outputsignal = Average_Nuc/Average_Cyto;
    case 2
        switch signaltype
            case 1
                outputsignal = mean(mean(mini_ratioIm(nuc_Y,nuc_X)));
            case 2
                outputsignal = mean(mean(cell_ch1(nuc_Y,nuc_X)));
            case 3
                outputsignal = mean(mean(cell_ch2(nuc_Y,nuc_X)));
            case 4
                outputsignal = mean(mean(cell_ch3(nuc_Y,nuc_X)));
        end        
    case 3
        switch signaltype
            case 1
                outputsignal = mean(mean(mini_ratioIm(cyto_Y,cyto_X)));
            case 2
                outputsignal = mean(mean(cell_ch1(cyto_Y,cyto_X)));
            case 3
                outputsignal = mean(mean(cell_ch2(cyto_Y,cyto_X)));
            case 4
                outputsignal = mean(mean(cell_ch3(cyto_Y,cyto_X)));
        end        
    case 4
        switch signaltype
            case 1
                outputsignal = mean(mean(mini_ratioIm(cell_Y,cell_X)));
            case 2
                outputsignal = mean(mean(cell_ch1(cell_Y,cell_X)));
            case 3
                outputsignal = mean(mean(cell_ch2(cell_Y,cell_X)));
            case 4
                outputsignal = mean(mean(cell_ch3(cell_Y,cell_X)));
        end        
end

if isnan(outputsignal)
    outputsignal=0;
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
savefile = ['AllFRETsignals_r' num2str(row) 'c' num2str(col) 'f' num2str(field)];
legendList=cell(1);
if exist([savefile '.mat'],'file')
    load(savefile);
    scell = str2num(get(handles.edit_cellNo,'String'));
    PosTime = find(mysignal1{scell}(:,2));
    if ~isempty(mysignal1{scell})
        figure;
        signalNames = get(handles.popupmenu_regionVar1,'String');
        hold on;plot(mysignal1{scell}(PosTime,1)/60,mysignal1{scell}(PosTime,2)/median(mysignal1{scell}(PosTime,2)),'b');
        s1Loc = get(handles.popupmenu_regionVar1,'Value');
        legendList{1} = signalNames{s1Loc};
        if get(handles.checkbox_variable2,'Value')
            hold on;plot(mysignal2{scell}(PosTime,1)/60,mysignal2{scell}(PosTime,2)/median(mysignal2{scell}(PosTime,2)),'r');
            s2Loc = get(handles.popupmenu_regionVar2,'Value');
            legendList{2} = signalNames{s2Loc};
        end
        if get(handles.checkbox_variable3,'Value')
            hold on;plot(mysignal3{scell}(PosTime,1)/60,mysignal3{scell}(PosTime,2)/median(mysignal3{scell}(PosTime,2)),'g');
            s3Loc = get(handles.popupmenu_regionVar3,'Value');
            legendList{3} =  signalNames{s3Loc};
        end
        if get(handles.checkbox_variable4,'Value')
            hold on;plot(mysignal4{scell}(PosTime,1)/60,mysignal4{scell}(PosTime,2)/median(mysignal4{scell}(PosTime,2)),'k');
            s4Loc = get(handles.popupmenu_regionVar4,'Value');
            legendList{4} = signalNames{s4Loc};
        end
        legend(legendList);
        xlabel('Time(hour)');
    end
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
templateCH= handles.templateCH;
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

maskOUTfilename = ['maskOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(templateCH)];
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
        MaskGenerating(3,[templateCH nominCH denominCH 1 last_tp],...
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
[filename,PathName,FilterIndex] = uigetfile('*.nd', 'Choose metamorph ND file');
if FilterIndex~=0
    set(handles.edit_ndfilename,'String',filename);
    handles.ndfilename = filename;
    handles.ndpathname = PathName;
end
guidata(hObject, handles);  

function [notp stagePos stageName waveName] = readndfile(pathname,filename)
% Search for number of string matches per line.  
notp=-1;
stagePos = [];
stageName = [];
waveName = [];
currentF = pwd;
cd(pathname);

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
cd(currentF);
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
        
        tokens   = regexp(handles.stageName{get(hObject,'Value')}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
        if ~isempty(tokens)
            row = tokens{1}{1};
            col = tokens{1}{2};
            set(handles.edit_row,'String',row);
            set(handles.edit_col,'String',col);
        else
            set(handles.edit_row,'String',num2str(get(hObject,'Value')));
            set(handles.edit_col,'String','1');
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

channel = handles.templateCH;
filetype = handles.filetype;
blankformat = get(handles.edit_fileformatBG,'String');
row = str2num(get(handles.edit_row,'String'));
col = str2num(get(handles.edit_col,'String'));
field = str2num(get(handles.edit_field,'String'));
plane = str2num(get(handles.edit_plane,'String'));
channelnames = handles.channelnames;
totalCHs = str2num(get(handles.edit_totalCHs,'String'));
tp = str2num(get(handles.edit_currentFrame,'String'));
template = [];
switch filetype
    case 1
        
        filename = sprintf(blankformat,channel,tp);
        if exist(filename,'file')
            template = im2double(imread(filename));
        else
            set(handles.edit_commu,'String',[filename ' does not exist.']);
        end
    case 2
        if exist(signalformat,'file')
            template = im2double(imread(blankformat,'Index',totalCHs*(tp-1)+channel));
        else
            set(handles.edit_commu,'String',[signalformat ' does not exist.']);
        end
    case 3
        filename = sprintf(blankformat,channelnames{channel},tp);
        if exist(filename,'file');
            template = im2double(imread(filename));
        else
            set(handles.edit_commu,'String',[filename ' does not exist.']);
        end
end

if ~isempty(template)
    switch get(hObject,'Value')
        case get(hObject,'Max')
            handles.useblank = 1;
            load fftexecutiontimes;
            h = 1/273*[1 4 7 4 1;4 16 26 26 4;7 26 41 26 7;4 16 26 16 4;1 4 7 4 1];
            handles.gaussian = h;
            handles.smooth_opt = detbestlength2(FFTrv,FFTiv,IFFTiv,size(template),size(h),isreal(template),isreal(h));
        case get(hObject,'Min')
            handles.useblank = 0;
            handles.smooth_opt = [];
    end
    set(handles.edit_commu,'String',['Assigned ' filename ' for blank']);
    guidata(hObject, handles);
end

else
    set(handles.edit_commu,'String','Blank removed');
end


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
tp = get(handles.edit_currentFrame,'String');
channel = handles.templateCH;
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
tp = get(handles.edit_currentFrame,'String');
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
tp = get(handles.edit_currentFrame,'String');
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
tp = get(handles.edit_currentFrame,'String');
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
