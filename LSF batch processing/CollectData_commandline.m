function CollectData_commandline(filetype,targetfolder,row,col,field,tps,signalformat,blankformat,channelnames,totalCHs,templateCH)
load trackingparameters

warning('off', 'all');

currentF = pwd;
cd(targetfolder);


first_tp = tps(1);
last_tp = tps(end);
load fftexecutiontimes
h_gaussian = 1/273*[1 4 7 4 1;4 16 26 26 4;7 26 41 26 7;4 16 26 16 4;1 4 7 4 1];
smooth_opt = detbestlength2(FFTrv,FFTiv,IFFTiv,[image_height image_width],size(h_gaussian),1,1);
H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
fileattrib(H5filename,'+w');
cellpath_name = ['/field' num2str(field) '/cellpath'];
sisterList_name = ['/field' num2str(field) '/sisterList'];
bg_name = ['/field' num2str(field) '/bg'];
maskdatasetname = ['/field' num2str(field) '/segmentsCH' num2str(templateCH)];
signal_name = ['/field' num2str(field) '/' outputname];
timestamp_name = ['/field' num2str(field) '/timestamp'];
selectedcells_name = ['/field' num2str(field) '/selectedcells'];

cellpathinfo = h5info(H5filename, cellpath_name);
sisterListinfo = h5info(H5filename, sisterList_name);
bginfo = h5info(H5filename, bg_name);
maskinfo = h5info(H5filename, maskdatasetname);

fid = H5F.open(H5filename,'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
    H5F.close(fid);
    cellpath_mat = h5read(H5filename,cellpath_name,[1 1 1], [cellpathinfo.Dataspace.Size(1) cellpathinfo.Dataspace.Size(2) cellpathinfo.Dataspace.Size(3)]);

    for tp=first_tp:last_tp
        cellpath{tp} = cellpath_mat(:,:,tp);
    end
else
    return
end

fid = H5F.open(H5filename,'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,sisterList_name,'H5P_DEFAULT')
    H5F.close(fid);
    sisterList_mat = h5read(H5filename,sisterList_name,[1 1 1], [sisterListinfo.Dataspace.Size(1) sisterListinfo.Dataspace.Size(2) sisterListinfo.Dataspace.Size(3)]);

    for tp=first_tp:last_tp
        sisterList{tp} = sisterList_mat(:,:,tp);
    end
else
    return
end

fid = H5F.open(H5filename,'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,bg_name,'H5P_DEFAULT')
    H5F.close(fid);
    bg_mat = h5read(H5filename,bg_name,[1 1 1], [bginfo.Dataspace.Size(1) bginfo.Dataspace.Size(2) bginfo.Dataspace.Size(3)]);

    for tp=first_tp:last_tp
        bg{tp} = bg_mat(:,:,tp);
    end
end

% Initializing storage variables
mysignal = zeros(size(cellpath{last_tp},1),last_tp-first_tp+1,4);
selected_cells = h5read(H5filename,selectedcells_name);
timestamp = 1:(last_tp-first_tp+1);

for tp=first_tp:last_tp
    % Determine image capture time based on the template channel
    if filetype == 3
        filename = sprintf(signalformat,channelnames{templateCH},first_tp);
        first_info = imfinfo(filename);
        filename = sprintf(signalformat,channelnames{templateCH},tp);
        current_info = imfinfo(filename);
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

    clc;display(['Processing time point:' num2str(tp) ' of ' num2str(last_tp)]);

    
    % Determine signals
    CH1im = loadsignal(filetype,signalformat,blankformat,channelnames,totalCHs,CH1,tp,h_gaussian,smooth_opt);
    CH2im = loadsignal(filetype,signalformat,blankformat,channelnames,totalCHs,CH2,tp,h_gaussian,smooth_opt);
    CH3im = loadsignal(filetype,signalformat,blankformat,channelnames,totalCHs,CH3,tp,h_gaussian,smooth_opt);
    ratioIm = calculateFRET(filetype,signalformat,blankformat,channelnames,totalCHs,tp,nominCH,denominCH,bg,filterParam1,filterParam2,bgsize,signalShiftN,signalShiftD);
    
    imwidth = size(CH1im,2);
    imheight = size(CH1im,1);
    
    startind = double([1 tp 1 1 1]);
    countind = [maskinfo.Dataspace.Size(1) 1 3 maskinfo.Dataspace.Size(4) maskinfo.Dataspace.Size(5)];
    allmasks = permute(h5read(H5filename,maskdatasetname,startind, countind),[4 5 1 3 2]);
    
    
    nucmask  = double(allmasks(:,:,:,1));
    cellmask  = double(allmasks(:,:,:,2));
    cytomask  = double(allmasks(:,:,:,3));
    
    for cellNo=selected_cells'

        ind_cellpath = pos_path(cellpath,sisterList,cellNo,first_tp,last_tp,imheight,imwidth);
        
        if tp <= size(ind_cellpath,1)
            
            if ~isempty(find(nucmask(:,:,cellNo),1))
                
                if cellsize~=(size(nucmask(:,:,cellNo),1)-1)/2
                    cellsize = (size(nucmask(:,:,cellNo),1)-1)/2;
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
                
                NucMask  = nucmask(:,:,cellNo) ;
                
                if ~isempty(find(cellmask(:,:,cellNo),1))
                    CellMask = cellmask(:,:,cellNo);
                else
                    CellMask = bwmorph(NucMask,'dilate',cytosize);
                    cellmask(:,:,cellNo) = CellMask;
                end
                
                if ~isempty(find(cytomask(:,:,cellNo),1))
                    CytoMask = cytomask(:,:,cellNo);
                else
                    CytoMask = CellMask-NucMask;
                    cytomask(:,:,cellNo) = CytoMask;
                end
                
                %for Nuclei region
                [nuc_Y, nuc_X] = find(NucMask);
                %for Cytosol region
                [cyto_Y, cyto_X] = find(CytoMask);
                %for Cell region
                [cell_Y, cell_X] = find(CellMask);
                
                clear NucMask CytoMask CellMask; 
                
                mysignal(cellNo,tp-first_tp+1,1) = signalOutputing(regiontype(1),signaltype(1),mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y);
                mysignal(cellNo,tp-first_tp+1,2) = signalOutputing(regiontype(2),signaltype(2),mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y);
                mysignal(cellNo,tp-first_tp+1,3) = signalOutputing(regiontype(3),signaltype(3),mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y);
                mysignal(cellNo,tp-first_tp+1,4) = signalOutputing(regiontype(4),signaltype(4),mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y);
                
                clear nuc_X nuc_Y cyto_X cyto_Y cell_X cell_Y cell_template mini_ratioIm cell_ch1 cell_ch2 cell_ch3
            end
        end
    end
    
    timestamp(tp-first_tp+1) = timestep;
    
    clear  ratioIm CH1im CH2im CH3im
end

fileattrib(H5filename,'+w');
% 
fid = H5F.open(H5filename,'H5F_ACC_RDWR','H5P_DEFAULT');

if ~H5L.exists(fid,signal_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' signal_name]);
else
    H5L.delete(fid,signal_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' signal_name]);
    H5F.close(fid);
end

h5create(H5filename, signal_name, [size(mysignal,1), size(mysignal,2), 4], 'Datatype', 'double', 'ChunkSize', [size(mysignal,1), size(mysignal,2), 1], 'Deflate', 9);
h5write(H5filename, signal_name, mysignal, [1 1 1], [size(mysignal,1) size(mysignal,2) 4]);

fid = H5F.open(H5filename,'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,timestamp_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' timestamp_name]);
else
    H5L.delete(fid,timestamp_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' timestamp_name]);
    H5F.close(fid);
end

h5create(H5filename, timestamp_name, last_tp-first_tp+1, 'Datatype', 'double');
h5write(H5filename, timestamp_name, timestamp);
cd(currentF);

display(['Successfully collected signals from frame' num2str(tps(1)) ':' num2str(tps(end)) ' and saved in ' H5filename]);


function outputsignal = signalOutputing(regiontype,signaltype,mini_ratioIm,cell_ch1,cell_ch2,cell_ch3,nuc_X,nuc_Y,cyto_X,cyto_Y,cell_X,cell_Y)
switch regiontype
    case 1
        switch signaltype
            case 4
                Average_Nuc = mean(mean(mini_ratioIm(nuc_Y,nuc_X)));
                Average_Cyto = mean(mean(mini_ratioIm(cyto_Y,cyto_X)));

            case 1
                Average_Nuc = mean(mean(cell_ch1(nuc_Y,nuc_X)));
                Average_Cyto = mean(mean(cell_ch1(cyto_Y,cyto_X)));
                
            case 2
                Average_Nuc = mean(mean(cell_ch2(nuc_Y,nuc_X)));
                Average_Cyto = mean(mean(cell_ch2(cyto_Y,cyto_X)));
            case 3
                Average_Nuc = mean(mean(cell_ch3(nuc_Y,nuc_X)));
                Average_Cyto = mean(mean(cell_ch3(cyto_Y,cyto_X)));
        end
        outputsignal = Average_Nuc/Average_Cyto;
    case 2
        switch signaltype
            case 4
                outputsignal = mean(mean(mini_ratioIm(nuc_Y,nuc_X)));
            case 1
                outputsignal = mean(mean(cell_ch1(nuc_Y,nuc_X)));
            case 2
                outputsignal = mean(mean(cell_ch2(nuc_Y,nuc_X)));
            case 3
                outputsignal = mean(mean(cell_ch3(nuc_Y,nuc_X)));
        end        
    case 3
        switch signaltype
            case 4
                outputsignal = mean(mean(mini_ratioIm(cyto_Y,cyto_X)));
            case 1
                outputsignal = mean(mean(cell_ch1(cyto_Y,cyto_X)));
            case 2
                outputsignal = mean(mean(cell_ch2(cyto_Y,cyto_X)));
            case 3
                outputsignal = mean(mean(cell_ch3(cyto_Y,cyto_X)));
        end        
    case 4
        switch signaltype
            case 4
                outputsignal = mean(mean(mini_ratioIm(cell_Y,cell_X)));
            case 1
                outputsignal = mean(mean(cell_ch1(cell_Y,cell_X)));
            case 2
                outputsignal = mean(mean(cell_ch2(cell_Y,cell_X)));
            case 3
                outputsignal = mean(mean(cell_ch3(cell_Y,cell_X)));
        end        
end

if isnan(outputsignal)
    outputsignal=0;
end

function ratioIm = calculateFRET(filetype,signalformat,blankformat,channelnames,totalCHs,tp,nominCH,denominCH,bg,filterParam1,filterParam2,bgsize,signalShiftN,signalShiftD)

nominIM = loadsignal(filetype,signalformat,blankformat,channelnames,totalCHs,nominCH,tp,h_gaussian,smooth_opt);
if denominCH == -1
    denomIM = im2double(ones(size(nominIM)));
else
    denomIM = loadsignal(filetype,signalformat,blankformat,channelnames,totalCHs,denominCH,tp,h_gaussian,smooth_opt);
end

% Load blank images if user chose to use BG points from BLANK images
if bgnomin_LOG==2
    nominBLK = loadblank(filetype,blankformat,channelnames,totalCHs,nominCH,tp);
    
end
if bgdenomin_LOG==2
    if denominCH == -1
        denomBLK = im2double(ones(size(nominIM)));
    else
        denomBLK = loadblank(filetype,blankformat,channelnames,totalCHs,denominCH,tp);
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

ratioIm =  normN./normD;


function outputim = loadsignal(filetype,signalformat,blankformat,channelnames,totalCHs,channel,tp,h_gaussian,smooth_opt)
outputim = [];



if useblank_LOG
    
    switch filetype
        case 1
            
            signal_filename = sprintf(signalformat,channel,tp);
            blank_filename = sprintf(blankformat,channel,tp);
            
            if exist(signal_filename,'file') 
                signalim = im2double(imread(signal_filename));
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
            if exist(signalformat,'file')
                outputim = im2double(imread(signalformat,'Index',totalCHs*(tp-1)+channel));
            else
                outputim = [];
            end
        case 3
            signal_filename = sprintf(signalformat,channelnames{channel},tp);
            blank_filename = sprintf(blankformat,channelnames{channel},tp);
            if exist(signal_filename,'file') 
                signalim = im2double(imread(signal_filename));
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
            
    end
else
    switch filetype
        case 1
            
            filename = sprintf(signalformat,channel,tp);
            if exist(filename,'file')
                outputim = im2double(imread(filename));
            else
                outputim = [];
            end
        case 2
            if exist(signalformat,'file')
                outputim = im2double(imread(signalformat,'Index',totalCHs*(tp-1)+channel));
            else
                outputim = [];
            end
        case 3
            filename = sprintf(signalformat,channelnames{channel},tp);
            if exist(filename,'file');
                outputim = im2double(imread(filename));
            else
                outputim = [];
            end
    end
end

function outputim = loadblank(filetype,blankformat,channelnames,totalCHs,channel,tp)

outputim = [];

switch filetype
    case 1
        
        filename = sprintf(blankformat,channel,tp);
        if exist(filename,'file')
            outputim = im2double(imread(filename));
        else
            outputim = [];
        end
    case 2
        if exist(signalformat,'file')
            outputim = im2double(imread(blankformat,'Index',totalCHs*(tp-1)+channel));
        else
            outputim = [];
        end
    case 3
        filename = sprintf(blankformat,channelnames{channel},tp);
        if exist(filename,'file');
            outputim = im2double(imread(filename));
        else
            outputim = [];
        end
end
