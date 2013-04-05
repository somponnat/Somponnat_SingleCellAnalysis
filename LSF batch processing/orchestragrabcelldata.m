function orchestragrabcelldata()

rows = 2;
cols = 4;
fields = 1:4;
channels = 4;
plane = 1;
nominCH = 2;
denominCH = 1;
tps = [1 130];
lind = 1;

matlabpool open;

for row = rows
    for col = cols
        for field = fields
            for templateCH = channels
                loopparam(lind,:) = [row col field templateCH];
                lind = lind+1;
            end
        end
    end
end

parfor lind = 1:size(loopparam,1)
    row = loopparam(lind,1);
    col = loopparam(lind,2);
    field = loopparam(lind,3);
    templateCH = loopparam(lind,4);
    collectData(row,col,field,plane,templateCH,nominCH,denominCH,tps);
end

matlabpool close;

function collectData(row,col,field,plane,templateCH,nominCH,denominCH,tps)
% hObject    handle to pushbutton_collect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

celltrackfile = ['celltrackOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(templateCH)];

% check if necessary file exist, if not quit jobs
if ~exist([celltrackfile '.mat'],'file')
    display([celltrackfile '.mat does not exist.']);
    return
end
% file exist, load the file
load(celltrackfile,'cellpath','sisterList','bg','mark');
if ~exist('cellpath','var') | ~exist('bg','var') | ~exist('sisterList','var')
    display([savefile '.mat does not contain necessary variables.']);
    return
end
%------------------------------
% Adjustable parameters
tpframe = 5; %minutes
filterParam1 = 2;
filterParam2 = 2;
displaygateL = 0.75;
displaygateH = 1.4;
signalShiftN = 2^8;
signalShiftD = 2^8;
cellsize = 80; %must be even integer
bgsize = 30;    %must be even integer
illumlogic = 0;
bglogic = 1;
maxI = 1e11;
invertLog = 0;
def_xshift = 0;
def_yshift = 0;
def_squaresize = 4;
def_cytosize = 4;
collectiontypeind = 3;
%------------------------------
% Edit this to assign which cells to collect data
collectedCells = 1:size(cellpath{tps(end)},1);
%------------------------------

% Initializing segmentation location info.
cytosize = cell(size(cellpath{tps(end)},1),1);
squaresize = cell(size(cellpath{tps(end)},1),1);
collectiontype = cell(size(cellpath{tps(end)},1),1);
shifting = cell(size(cellpath{tps(end)},1),1);

segmentfile = ['datacollection_r' num2str(row) 'c' num2str(col) 'f' num2str(field)];

if exist([segmentfile '.mat'],'file')
    load(segmentfile,'shifting','squaresize','cytosize','collectiontype');
end


%------------------------------
fileformat = 'r%02.0fc%02.0ff%02.0fp%02.0frc%1.0f-ch1sk%ufk1fl1.tiff';
fftw('planner', 'hybrid');
warning('off', 'all');

filename = sprintf(fileformat,row,col,field,plane,1,tps(1));
temp = im2double(imread(filename));

for tp=tps(1):tps(end)
    clc;display([celltrackfile ' - Currently processing frame: ' num2str(tp) ' of ' num2str(tps(end)-tps(1)+1)]);
    
    filename = sprintf(fileformat,row,col,field,plane,templateCH,tp);
    template = imread(filename);
    
    if nominCH == -1
        nominIM = double(ones(size(temp)));
    else
        filename = sprintf(fileformat,row,col,field,plane,nominCH,tp);
        nominIM = double(imread(filename));
    end
    
    if denominCH == -1
        denomIM = double(ones(size(temp)));
    else
        filename = sprintf(fileformat,row,col,field,plane,denominCH,tp);
        denomIM = double(imread(filename));
    end
    
    
 
    if illumlogic
        normN = double(ifft2(ifftshift(fftshift(fft2(nominIM)).*hbutter(nominIM,filterParam1,filterParam2))));
        normD = double(ifft2(ifftshift(fftshift(fft2(denomIM)).*hbutter(denomIM,filterParam1,filterParam2))));
    else
        normN = nominIM;
        normD = denomIM;
    end
    clear nominIM denomIM;
    templatesize = size(template);
    if ~isempty(bg)
        for b=1:size(bg{tp},1)
            xL=max(bg{tp}(b,1)-bgsize/2,1);
            xR=min(bg{tp}(b,1)+bgsize/2,templatesize(2));
            yL=max(bg{tp}(b,2)-bgsize/2,1);
            yR=min(bg{tp}(b,2)+bgsize/2,templatesize(1));
            BG_N(b) = mean(mean(normN(yL:yR,xL:xR)));
            BG_D(b) = mean(mean(normD(yL:yR,xL:xR)));
        end
        meanBG_normN = mean(BG_N);
        meanBG_normD = mean(BG_D);
        
        if bglogic
            normN = normN-meanBG_normN;
            normD = normD-meanBG_normD;
       end
        
    end

    normN = normN+signalShiftN;
    normD = normD+signalShiftD;
    
%     ratioIm =  mat2gray(normN./normD,[displaygateL displaygateH]);
%     figure(1),imshow(ratioIm,[]); hold on;
    %plot(cellpath{tp}(:,1),cellpath{tp}(:,2),'xr');
    for ccell=1:size(cellpath{tp},1)
        cellInd = find(ccell==collectedCells);
        if ~isempty(cellInd)
            xL=max(cellpath{tp}(ccell,1)-cellsize/2,1);
            xR=min(cellpath{tp}(ccell,1)+cellsize/2-1,templatesize(2));
            yL=max(cellpath{tp}(ccell,2)-cellsize/2,1);
            yR=min(cellpath{tp}(ccell,2)+cellsize/2-1,templatesize(1));
            
            cell_template = template(yL:yR,xL:xR);
            mini_normN    = normN(yL:yR,xL:xR);
            mini_normD    = normD(yL:yR,xL:xR);
            
            % determine pixel coordinate for different cell regions
            if isempty(cytosize{ccell})
                mycytosize = def_cytosize;
            else
                mycytosize = cytosize{ccell};
            end
            [x1 y1 nucMask cytoMask] = templateToCentroid(cell_template,cellsize/2,cellsize/2,mycytosize,maxI,invertLog);
            size_cell_template = size(cell_template);
            
            %for Nuclei region
            [nuc_Y nuc_X] = find(nucMask);
            
            %for Cytosol region
            [cyto_Y cyto_X] = find(cytoMask);
            
            %for Box region            
            sxL=max(def_xshift+cellsize/2-def_squaresize,1);
            sxR=min(def_xshift+cellsize/2+def_squaresize-1,size_cell_template(2));
            syL=max(def_yshift+cellsize/2-def_squaresize,1);
            syR=min(def_yshift+cellsize/2+def_squaresize-1,size_cell_template(1));
            
            [Xv Yv] = meshgrid(sxL:sxR,syL:syR);
            box_X = Xv(:);
            box_Y = Yv(:);
            
            clear Xv Yv sxL sxR syL syR;
            % for multi
            if ~isempty(collectiontype{ccell}) & ~isempty(squaresize{ccell}) & ~isempty(shifting{ccell})
                switch collectiontype{ccell}
                    case 1
                        multi_X = nuc_X;
                        multi_Y = nuc_Y;

                    case 2
                        multi_X = cyto_X;
                        multi_Y = cyto_Y;

                    case 3
                        def_squaresize = squaresize{ccell};
                        def_xshift = shifting{ccell}(1);
                        def_yshift = shifting{ccell}(2);
                        
                        sxL=max(def_xshift+cellsize/2-def_squaresize,1);
                        sxR=min(def_xshift+cellsize/2+def_squaresize-1,size_cell_template(2));
                        syL=max(def_yshift+cellsize/2-def_squaresize,1);
                        syR=min(def_yshift+cellsize/2+def_squaresize-1,size_cell_template(1));
                        
                        [Xv Yv] = meshgrid(sxL:sxR,syL:syR);
                        multi_X = Xv(:);
                        multi_Y = Yv(:);
                        clear Xv Yv sxL sxR syL syR;
                end
            else
                switch collectiontypeind
                    case 1
                        multi_X = nuc_X;
                        multi_Y = nuc_Y;
                    case 2
                        multi_X = cyto_X;
                        multi_Y = cyto_Y;
                    case 3
                        
                        multi_X = box_X;
                        multi_Y = box_Y;
                        
                end
            end
            mysignal1(tp,ccell) = signalOutputing(1,mini_normN,mini_normD,multi_X,multi_Y,nuc_X,nuc_Y,cyto_X,cyto_Y,box_X,box_Y);
            mysignal2(tp,ccell) = signalOutputing(2,mini_normN,mini_normD,multi_X,multi_Y,nuc_X,nuc_Y,cyto_X,cyto_Y,box_X,box_Y);
            mysignal3(tp,ccell) = signalOutputing(3,mini_normN,mini_normD,multi_X,multi_Y,nuc_X,nuc_Y,cyto_X,cyto_Y,box_X,box_Y);
            mysignal4(tp,ccell) = signalOutputing(4,mini_normN,mini_normD,multi_X,multi_Y,nuc_X,nuc_Y,cyto_X,cyto_Y,box_X,box_Y);
        end
    end
    clear normD normN denomIM nominIM mini_normN mini_normD  
    clc;
    
end

savefile = ['AllFRETsignals_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_ch' num2str(templateCH)];
timestep = tpframe*((tps(1):tps(end))-1); %minutes
save(savefile,'timestep','-v7.3');
save(savefile,'mysignal1','-v7.3','-append');
save(savefile,'mysignal2','-v7.3','-append');
save(savefile,'mysignal3','-v7.3','-append');
save(savefile,'mysignal4','-v7.3','-append');
display(['done collecting and saving data to ' savefile '.mat']);

function outputsignal = signalOutputing(regiontype,mini_normN,mini_normD,multi_X,multi_Y,nuc_X,nuc_Y,cyto_X,cyto_Y,box_X,box_Y)

switch regiontype
    case 1
        
        Average_N = mean(mean(mini_normN(multi_Y,multi_X)));
        Average_D = mean(mean(mini_normD(multi_Y,multi_X)));
        outputsignal = Average_N/Average_D;
        
    case 2
        
        Average_N = mean(mean(mini_normN(nuc_Y,nuc_X)));
        Average_D = mean(mean(mini_normD(nuc_Y,nuc_X)));
        outputsignal = Average_N/Average_D;
        
    case 3
        
        Average_N = mean(mean(mini_normN(cyto_Y,cyto_X)));
        Average_D = mean(mean(mini_normD(cyto_Y,cyto_X)));
        outputsignal = Average_N/Average_D;
        
    case 4
        
        Average_N = mean(mean(mini_normN(box_Y,box_X)));
        Average_D = mean(mean(mini_normD(box_Y,box_X)));
        outputsignal = Average_N/Average_D;
        
end



function [x y nucMask cytoMask] = templateToCentroid(M,xg,yg,mycytosize,maxI,invertLog)
M  = imadjust(M);
if invertLog
    inverted = ((maxI)-1)-M;
end
BWc = zeros(size(M));
for i=2:12
    edgedIm = edge(M,'canny',0,i);
    BW = imfill(edgedIm,'holes');
    
    BW = bwmorph(BW,'open',1);
    BW = bwselect(BW,xg,yg);
    
    BWc = BWc | BW;

end

BW = BWc;
S  = regionprops(BW, 'centroid');

if isempty(find(BW==0)) | isempty(find(BW==1))
    x = xg;
    y = yg;

else
    x = round(S.Centroid(1));
    y = round(S.Centroid(2));
end

nucMask = bwselect(BW,xg,yg);
BW3 = bwmorph(nucMask,'dilate',mycytosize);
cytoMask = BW3-nucMask;


function out = hbutter(im,d,n)
height = size(im,1);
width = size(im,2);
[x,y] = meshgrid(-floor(width/2):floor((width-1)/2),-floor(height/2):floor((height-1)/2));
out = 1-(1./(1+(sqrt(x.^2+y.^2)/d).^(2*n)));