function CollectData_commandline(filetype,ndpathname,row,col,field,tps,signalformat,blankformat,channelnames,totalCHs)
load trackingparameters

currentPath = pwd;
eval('cd ..');
addpath(genpath([pwd filesep 'ThirdParty']),'-end');
cd(currentPath);

first_tp = tps(1);
last_tp = tps(end);
load fftexecutiontimes
h_gaussian = 1/273*[1 4 7 4 1;4 16 26 26 4;7 26 41 26 7;4 16 26 16 4;1 4 7 4 1];
smooth_opt = detbestlength2(FFTrv,FFTiv,IFFTiv,[image_height image_width],size(h_gaussian),1,1);
H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
fileattrib(fullfile(ndpathname,H5filename),'+w');
cellpath_name = ['/field' num2str(field) '/cellpath'];
selectedcells_name = ['/field' num2str(field) '/selectedcells'];
signal_name = ['/field' num2str(field) '/outputsignal' num2str(outputsignalNo)];
timestamp_name = ['/field' num2str(field) '/timestamp' num2str(outputsignalNo)];

cellpathinfo = h5info(fullfile(ndpathname,H5filename), cellpath_name);

fileattrib(fullfile(ndpathname,H5filename),'+w');

fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
    H5F.close(fid);
    cellpath_mat = h5read(fullfile(ndpathname,H5filename),cellpath_name,[1 1 1], [cellpathinfo.Dataspace.Size(1) cellpathinfo.Dataspace.Size(2) cellpathinfo.Dataspace.Size(3)]);
else
    return
end

basename='outputsignal';
info = h5info(fullfile(ndpathname,H5filename),['/field' num2str(field)]);
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
mysignal = zeros(size(cellpath_mat,3),last_tp-first_tp+1,4);
h5create(fullfile(ndpathname,H5filename), signal_name, [size(mysignal,1), size(mysignal,2), 4], 'Datatype', 'double', 'ChunkSize', [size(mysignal,1), size(mysignal,2), 1], 'Deflate', 9);
h5write(fullfile(ndpathname,H5filename), signal_name, mysignal, [1 1 1], [size(mysignal,1) size(mysignal,2) 4]);
h5writeatt(fullfile(ndpathname,H5filename),signal_name,'outputsignal_name',output_name);

% mysignal1 labeling
region1 = regions{regiontype(1)};
signal1 = signals{signaltype(1)};
h5writeatt(fullfile(ndpathname,H5filename),signal_name,'signal1',[region1 '-' signal1]);
% mysignal2 labeling
region2 = regions{regiontype(2)};
signal2 = signals{signaltype(2)};
h5writeatt(fullfile(ndpathname,H5filename),signal_name,'signal2',[region2 '-' signal2]);
% mysignal3 labeling
region3 = regions{regiontype(3)};
signal3 = signals{signaltype(3)};
h5writeatt(fullfile(ndpathname,H5filename),signal_name,'signal3',[region3 '-' signal3]);
% mysignal4 labeling
region4 = regions{regiontype(4)};
signal4 = signals{signaltype(4)};
h5writeatt(fullfile(ndpathname,H5filename),signal_name,'signal4',[region4 '-' signal4]);


timestamp = first_tp:(last_tp-first_tp+1);
fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,timestamp_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' timestamp_name]);
else
    H5L.delete(fid,timestamp_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' timestamp_name]);
    H5F.close(fid);
end

h5create(fullfile(ndpathname,H5filename), timestamp_name, last_tp-first_tp+1, 'Datatype', 'double');
h5write(fullfile(ndpathname,H5filename), timestamp_name, timestamp);

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
sSerpen_idx = find(serpen(:,1)> size(cellpath_mat,3)*0.5);
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


for tp=first_tp:last_tp
    % Determine image capture time based on the template channel
    if filetype == 3
        filename = sprintf(signalformat,channelnames{cellCH},first_tp);
        first_info = imfinfo(fullfile(ndpathname,filename));
        filename = sprintf(signalformat,channelnames{cellCH},tp);
        current_info = imfinfo(fullfile(ndpathname,filename));
    end

    switch filetype
        case 3
            [~, ~, D, H, MN, S]  = datevec(datenum(current_info.DateTime,'yyyymmdd HH:MM:SS.FFF')-datenum(first_info.DateTime,'yyyymmdd HH:MM:SS.FFF'));
            hour = 24*D+H;
            minute = MN;
            second = round(S);
        otherwise
            hour = floor(1/60*(tp-1)*( minstamp + secstamp/60 )); %#ok<*ST2NM>
            minute = mod(floor((tp-1)*( minstamp + secstamp/60 )),60);
            second = mod((tp-1)*( minstamp*60 + secstamp),60);
    end
    
    timestep = hour*60+ minute + second/60;
    h5write(fullfile(ndpathname,H5filename), timestamp_name, timestep,tp-first_tp+1,1);
    
    display(['Processing time point:' num2str(tp) ' of ' num2str(last_tp) ',' H5filename]);

    
    % Determine signals
    CH1im = loadsignal(filetype,ndpathname,signalformat,blankformat,channelnames,totalCHs,CH1,tp,h_gaussian,smooth_opt);
    CH2im = loadsignal(filetype,ndpathname,signalformat,blankformat,channelnames,totalCHs,CH2,tp,h_gaussian,smooth_opt);
    ratioim = calculateFRET(filetype,ndpathname,signalformat,blankformat,channelnames,totalCHs,tp,[],filterParam1,filterParam2,bgsize,signalShiftN,signalShiftD,h_gaussian,smooth_opt,useblank_LOG);

    nuc_im   = imread(fullfile(ndpathname,nucFolder, sprintf(signalformat,channelnames{nucCH},tp)));
    if exist(fullfile(ndpathname,cellFolder,sprintf(signalformat,channelnames{cellCH},tp)),'file')
        cell_im  = imread(fullfile(ndpathname,cellFolder,sprintf(signalformat,channelnames{cellCH},tp)));
    else
        cell_im  = imread(fullfile(ndpathname,cellFolder,sprintf(signalformat,channelnames{nucCH},tp)));
    end
    
    for s=sSerpen_idx'
        pos_time = (serpen(s,3):serpen(s,4));
        if ismember(tp,pos_time)
            cellNo = serpen(s,2);
            cX = cellpath_mat(cellNo,1,tp);
            cY = cellpath_mat(cellNo,2,tp);

            NucMask  = bwselect(nuc_im,cX,cY,8);
            CellMask = bwselect(cell_im,cX,cY,8);
            CytoMask = CellMask-NucMask;
            
            %for Nuclei region
            [nuc_Y, nuc_X] = find(NucMask);
            %for Cytosol region
            [cyto_Y, cyto_X] = find(CytoMask);
            %for Cell region
            [cell_Y, cell_X] = find(CellMask);
            
            mysignal(cellNo,tp-first_tp+1,1) = signalOutputing(regiontype(1),signaltype(1),CH1im,CH2im,ratioim,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y);
            mysignal(cellNo,tp-first_tp+1,2) = signalOutputing(regiontype(2),signaltype(2),CH1im,CH2im,ratioim,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y);
            mysignal(cellNo,tp-first_tp+1,3) = signalOutputing(regiontype(3),signaltype(3),CH1im,CH2im,ratioim,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y);
            mysignal(cellNo,tp-first_tp+1,4) = signalOutputing(regiontype(4),signaltype(4),CH1im,CH2im,ratioim,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y);
            clear nuc_X nuc_Y cyto_X cyto_Y cell_X cell_Y
        end
    end
    h5write(fullfile(ndpathname,H5filename), signal_name, mysignal(:,tp-first_tp+1,:), [1 tp-first_tp+1 1], [size(mysignal,1) 1 4]);
    clear  CH1im CH2im nuc_im cell_im ratioim
end

display(['Successfully collected signals from frame' num2str(tps(1)) ':' num2str(tps(end)) ' and saved in ' H5filename]);

function outputsignal = signalOutputing(regiontype,signaltype,CH1im,CH2im,ratioim,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y)
nucData = [];
cytoData = [];
cellData = [];
switch regiontype
    case 1
        switch signaltype

            case 1
                for i=1:length(nuc_Y)
                    nucData(i) = CH1im(nuc_Y(i),nuc_X(i));
                end

                for i=1:length(cyto_Y)
                    cytoData(i) = CH1im(cyto_Y(i),cyto_X(i));
                end
                
                if ~isempty(nucData)
                    Average_Nuc = median(nucData);
                else
                    Average_Nuc = 0;
                end
                if ~isempty(cytoData)
                    Average_Cyto = median(cytoData);
                else
                    Average_Cyto = 0;
                end
                
            case 2
                for i=1:length(nuc_Y)
                    nucData(i) = CH2im(nuc_Y(i),nuc_X(i));
                end

                for i=1:length(cyto_Y)
                    cytoData(i) = CH2im(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(nucData)
                    Average_Nuc = median(nucData);
                else
                    Average_Nuc = 0;
                end
                if ~isempty(cytoData)
                    Average_Cyto = median(cytoData);
                else
                    Average_Cyto = 0;
                end

            case 3
                for i=1:length(nuc_Y)
                    nucData(i) = ratioim(nuc_Y(i),nuc_X(i));
                end

                for i=1:length(cyto_Y)
                    cytoData(i) = ratioim(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(nucData)
                    Average_Nuc = median(nucData);
                else
                    Average_Nuc = 0;
                end
                if ~isempty(cytoData)
                    Average_Cyto = median(cytoData);
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
                    nucData(i) = CH1im(nuc_Y(i),nuc_X(i));
                end
                if ~isempty(nucData)
                    outputsignal = median(nucData);
                else
                    outputsignal = 0;
                end
            case 2
                for i=1:length(nuc_Y)
                    nucData(i) = CH2im(nuc_Y(i),nuc_X(i));
                end
                if ~isempty(nucData)
                    outputsignal = median(nucData);
                else
                    outputsignal = 0;
                end
            case 3
                for i=1:length(nuc_Y)
                    nucData(i) = ratioim(nuc_Y(i),nuc_X(i));
                end
                if ~isempty(nucData)
                    outputsignal = median(nucData);
                else
                    outputsignal = 0;
                end
        end        
    case 3
        switch signaltype
            case 1
                for i=1:length(cyto_Y)
                    cytoData(i) = CH1im(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(cytoData)
                    outputsignal = median(cytoData);
                else
                    outputsignal = 0;
                end
            case 2
                for i=1:length(cyto_Y)
                    cytoData(i) = CH2im(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(cytoData)
                    outputsignal = median(cytoData);
                else
                    outputsignal = 0;
                end
            case 3
                for i=1:length(cyto_Y)
                    cytoData(i) = ratioim(cyto_Y(i),cyto_X(i));
                end
                if ~isempty(cytoData)
                    outputsignal = median(cytoData);
                else
                    outputsignal = 0;
                end
        end        
    case 4
        switch signaltype
            case 1
                for i=1:length(cell_Y)
                    cellData(i) = CH1im(cell_Y(i),cell_X(i));
                end
                if ~isempty(cellData)
                    outputsignal = median(cellData);
                else
                    outputsignal = 0;
                end
            case 2
                for i=1:length(cell_Y)
                    cellData(i) = CH2im(cell_Y(i),cell_X(i));
                end
                if ~isempty(cellData)
                    outputsignal = median(cellData);
                else
                    outputsignal = 0;
                end
            case 3
                for i=1:length(cell_Y)
                    cellData(i) = ratioim(cell_Y(i),cell_X(i));
                end
                if ~isempty(cellData)
                    outputsignal = median(cellData);
                else
                    outputsignal = 0;
                end
        end        
end

function outputim = loadsignal(filetype,ndpathname,signalformat,blankformat,channelnames,totalCHs,channel,tp,h_gaussian,smooth_opt)
outputim = [];
load trackingparameters

if useblank_LOG
    
    switch filetype
        case 1
            
            signal_filename = sprintf(signalformat,channel,tp);
            blank_filename = sprintf(blankformat,channel,tp);
            
            if exist(fullfile(ndpathname,signal_filename),'file') 
                signalim = im2double(imread(fullfile(ndpathname,signal_filename)));
                if exist(blank_filename,'file')
                    blankim = im2double(imread(blank_filename));
                    
                    smoothed = fftolamopt2(blankim,h_gaussian,smooth_opt,'same');
                    normim =(smoothed./mean(smoothed(:)));
                    clear blankim smoothed;
                    outputim = signalim./normim;
                else 

                    outputim = signalim;
                end
            else
                outputim = [];
            end
            
            
        case 2
            if exist(fullfile(ndpathname,signalformat),'file')
                outputim = im2double(imread(fullfile(ndpathname,signalformat),'Index',totalCHs*(tp-1)+channel));
            else
                outputim = [];
            end
        case 3
            signal_filename = sprintf(signalformat,channelnames{channel},tp);
            blank_filename = sprintf(blankformat,channelnames{channel},tp);
            if exist(fullfile(ndpathname,signal_filename),'file') 
                signalim = im2double(imread(fullfile(ndpathname,signal_filename)));
                if exist(fullfile(ndpathname,blank_filename),'file')
                    blankim = im2double(imread(fullfile(ndpathname,blank_filename)));
                    smoothed = fftolamopt2(blankim,h_gaussian,smooth_opt,'same');
                    normim =(smoothed./mean(smoothed(:)));
                    clear blankim smoothed;
                    outputim = signalim./normim;
                else 
                    outputim = signalim;
                end
            else
                outputim = [];
            end
            
    end
else
    switch filetype
        case 1
            
            filename = sprintf(signalformat,channel,tp);
            if exist(fullfile(ndpathname,filename),'file')
                outputim = im2double(imread(fullfile(ndpathname,filename)));
            else
                outputim = [];
            end
        case 2
            if exist(fullfile(ndpathname,signalformat),'file')
                outputim = im2double(imread(fullfile(ndpathname,signalformat),'Index',totalCHs*(tp-1)+channel));
            else
                outputim = [];
            end
        case 3
            filename = sprintf(signalformat,channelnames{channel},tp);
            if exist(fullfile(ndpathname,filename),'file');
                outputim = im2double(imread(fullfile(ndpathname,filename)));
            else
                outputim = [];
            end
    end
end

function outputim = loadblank(filetype,ndpathname,blankformat,channelnames,totalCHs,channel,tp)

outputim = [];

switch filetype
    case 1
        
        filename = sprintf(blankformat,channel,tp);
        if exist(fullfile(ndpathname,filename),'file')
            outputim = im2double(imread(filename));
        else
            outputim = [];
        end
    case 2
        if exist(fullfile(ndpathname,blankformat),'file')
            outputim = im2double(imread(fullfile(ndpathname,blankformat),'Index',totalCHs*(tp-1)+channel));
        else
            outputim = [];
        end
    case 3
        filename = sprintf(blankformat,channelnames{channel},tp);
        if exist(fullfile(ndpathname,filename),'file');
            outputim = im2double(imread(fullfile(ndpathname,filename)));
        else
            outputim = [];
        end
end

function out = hbutter(im,d,n)
height = size(im,1);
width = size(im,2);
[x,y] = meshgrid(-floor(width/2):floor((width-1)/2),-floor(height/2):floor((height-1)/2));
out = 1-(1./(1+(sqrt(x.^2+y.^2)/d).^(2*n)));

function ratioim = calculateFRET(filetype,ndpathname,signalformat,blankformat,channelnames,totalCHs,tp,bg,filterParam1,filterParam2,bgsize,signalShiftN,signalShiftD,h_gaussian,smooth_opt,useblank_LOG)
load trackingparameters
nominIM = loadsignal(filetype,ndpathname,signalformat,blankformat,channelnames,totalCHs,nominCH,tp,h_gaussian,smooth_opt);
if denominCH == -1
    denomIM = im2double(ones(size(nominIM)));
else
    denomIM = loadsignal(filetype,ndpathname,signalformat,blankformat,channelnames,totalCHs,denominCH,tp,h_gaussian,smooth_opt);
end

% Load blank images if user chose to use BG points from BLANK images
if bgnomin_LOG==2
    nominBLK = loadblank(filetype,ndpathname,blankformat,channelnames,totalCHs,nominCH,tp);
    
end
if bgdenomin_LOG==2
    if denominCH == -1
        denomBLK = im2double(ones(size(nominIM)));
    else
        denomBLK = loadblank(filetype,ndpathname,blankformat,channelnames,totalCHs,denominCH,tp);
    end
end

if illumlogic
    normN = ifft2(ifftshift(fftshift(fft2(nominIM)).*hbutter(nominIM,filterParam1,filterParam2)));
    normD = ifft2(ifftshift(fftshift(fft2(denomIM)).*hbutter(denomIM,filterParam1,filterParam2)));
else
    normN = nominIM;
    normD = denomIM;
end



switch bgnomin_LOG
    case 1
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
    case 2
        BG_N = median(nominBLK(:));
    otherwise
        BG_N = 0;
end

switch bgdenomin_LOG
    case 1
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
    case 2
        BG_D = median(denomBLK(:));
    otherwise
        BG_D = 0;
end


normN = normN-mean(BG_N);
normD = normD-mean(BG_D);

normN = normN+signalShiftN;
normD = normD+signalShiftD;

ratioim =  normN./normD;
