function maskgeneratingND_local()

clear all;
% Define information about input images and necessary parameters-----------
ndfilename = '130722.nd';
templateCH = 2;
cellCH1 = 1;
cellsize = 80; % should be at least 30
dilateSize = 10;
sourcefolder = 'Q:\sorger\data\NIC\Bernhard\130722';
%sourcefolder = 'Q:\sorger\data\NIC\Pat\02-03-2013\';
%sourcefolder = '~/files/ImStor/sorger/data/NIC/Pat/02-03-2013';
%------------------------------------------------
prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(sourcefolder,ndfilename);

tps = [1 notp];
sites = [17 19 20 21 22 24 37];%1:length(stagePos);

isOpen = matlabpool('size') > 0
if ~isOpen
    matlabpool local 
end


parfor i=1:length(sites)
    warning off;
    site = sites(i);
    fileformat = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
    %tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)f(?<field>\d+)','tokens');
    tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
    
    if ~isempty(tokens)
        row = str2num(tokens{1}{1});
        col = str2num(tokens{1}{2});
    else
        row = site;
        col = 1;
    end
    
    %field = tokens{1}{3};
    field = 1;
    plane = 1;
    maskgen_commandline(3,sourcefolder,row,col,field,plane,templateCH,cellCH1,tps,fileformat,channelnames,cellsize,dilateSize);
end


% This is the main function used to determine the segmentation
function [nucMask,cellMask] = templateToCentroid(nucOri,cellOri,cellCoord,neighborCoord,dilateSize,borderX,borderY)
nucMask = zeros(size(nucOri));
cellMask = zeros(size(nucOri));

nucOri = nucOri(borderY,borderX);
cellOri = cellOri(borderY,borderX);

allCoord = [cellCoord;neighborCoord];

BWc = zeros(size(nucOri));

% Loop through a range of sigma to better detect cell edges
for i=1.2:1:6
    
    % reduce edged image down to closed object
    edgedIm = bwperim(edge(imadjust(nucOri),'canny',0,i));
    % fill closed object
    BW = imfill(edgedIm,'holes');
    BW = bwmorph(BW,'open',1);
    
    % loop through both center cells and neighbor cells
    for j=1:size(allCoord,1)
        
        % only process cells within the assigned cropped image
        if allCoord(j,1)<=size(BW,2) && allCoord(j,2)<=size(BW,1) && allCoord(j,1)>0 && allCoord(j,2)>0
            
            % identify closed object at the current coordinate
            BW_ind = bwselect(BW,allCoord(j,1),allCoord(j,2));
            
            % calculate perimeter, area and majorAxisLength 
            S  = regionprops(BW_ind,'Perimeter','Area','MajorAxisLength');
            
            % only update the mask if parameters stay within the specified
            % constain
            if ~isempty(S) && S.MajorAxisLength < 30 &&  S.Area > 10 && 4*pi*S.Area/S.Perimeter^2 >0.9
                BWc = BWc | BW_ind;
            else
                
                % if no nuclear mask is detected with template channel,
                % attempt to detect nuclear segment using whole cell
                % channel
                BW_ind = bwperim(edge(imadjust(cellOri),'canny',0,i));
                BW_ind = imfill(BW_ind,'holes');
                BW_ind = bwmorph(BW_ind,'open',1);
                BW_ind = bwselect(BW_ind,allCoord(j,1),allCoord(j,2));
                S  = regionprops(BW_ind,'Perimeter','Area','MajorAxisLength');
                if ~isempty(S) && S.Area > 10 && 4*pi*S.Area/S.Perimeter^2 >0.5
                    BWc = BWc | BW_ind;    
                else
                    % if still cannot detect nucleus, imply use extended
                    % disk with the radius of 10 pixel
                    BW_ind = zeros(size(BWc));
                    se = strel('disk',10);
                    BW_ind(allCoord(j,2),allCoord(j,1)) = 1;
                    BW_ind = imdilate(BW_ind,se);
                    BWc = BWc | BW_ind;
                end
            end
            clear BW_ind
        end
        clear S;
    end
end

% now, attempt to detect the cell boundary by extending from th nuclear
% segment
BW = bwselect(BWc,allCoord(:,1),allCoord(:,2));
nucMask_temp = bwselect(BW,allCoord(1,1),allCoord(1,2));

mydistBW = bwdist(BW);
watershedBW  = watershed(mydistBW);

watershedNo = unique(watershedBW);
ThresholdedOrigImage = zeros(size(BW));
se = strel('disk',5);

for i=2:length(watershedNo)
    cellNo = watershedNo(i);
    temp = zeros(size(BW));
    selectedTemp = (watershedBW==cellNo);
    [r, c] = find(selectedTemp & ~BW);
    for j =1:length(r)
        temp(r(j),c(j)) = cellOri(r(j),c(j));
    end
    clear r c;
    % determine local threshold based on the watershed area
    ind_ThresImage = im2bw(cellOri,2*graythresh(temp));
    ind_ThresImage = imfill(ind_ThresImage,'holes');
    ind_ThresImage = bwmorph(ind_ThresImage,'open','Inf');
    ThresholdedOrigImage = ThresholdedOrigImage | ind_ThresImage;

end


if isempty(find(ThresholdedOrigImage,1))
    ThresholdedOrigImage = im2bw(cellOri,graythresh(cellOri));
end

PrimaryLabelImage = bwlabel(BW);

% use cellprofiler IdentifySecondaryObject by propagation
[labels_out,d]=IdentifySecPropagateSubfunction(PrimaryLabelImage,cellOri,ThresholdedOrigImage|BW,0.1);

secondaryLabel = imfill(labels_out);
selectedCell_group = secondaryLabel(allCoord(1,2),allCoord(1,1));
[Y,X] = find(secondaryLabel== selectedCell_group);
cellMask_temp = zeros(size(nucMask_temp));

% ensure that the cytosolic mask is at least 3 pixel in size 
if ~isempty(find(nucMask_temp,1))
    se = strel('disk',dilateSize);
    se2 = strel('disk',3);
    cellMask_Largest = imdilate(nucMask_temp,se);
    cellMask_smallest = imdilate(nucMask_temp,se2);
    for i=1:length(Y)
        cellMask_temp(Y(i),X(i)) = 1;
    end
    
    cellMask_temp = (cellMask_temp | cellMask_smallest) & cellMask_Largest ;
    
    cellMask_temp = bwselect(cellMask_temp,allCoord(1,1),allCoord(1,2));
end


nucMask(borderY,borderX) = nucMask_temp;
cellMask(borderY,borderX) = cellMask_temp;


function maskgen_commandline(filetype,targetfolder,row,col,field,plane,templateCH,cellCH1,tps,fileformat,channelnames,cellsize,dilateSize)

firsttp = tps(1);
lasttp = tps(end);

% This section specific the naming of HDF5 file and subfolder
H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
cellpath_name = ['/field' num2str(field) '/cellpath'];
sisterList_name = ['/field' num2str(field) '/sisterList'];
bg_name = ['/field' num2str(field) '/bg'];
imagesegment = ['/field' num2str(field) '/segmentsCH' num2str(templateCH)];
selectedcells_name = ['/field' num2str(field) '/selectedcells'];

cellpathinfo = h5info(fullfile(targetfolder,H5filename), cellpath_name);
sisterListinfo = h5info(fullfile(targetfolder,H5filename), sisterList_name);
bginfo = h5info(fullfile(targetfolder,H5filename), bg_name);

fileattrib(fullfile(targetfolder,H5filename),'+w');
fid = H5F.open(fullfile(targetfolder,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
    H5F.close(fid);
    
    % Read in the cellpath matrix
    cellpath_mat = h5read(fullfile(targetfolder,H5filename),cellpath_name,[1 1 1], [cellpathinfo.Dataspace.Size(1) cellpathinfo.Dataspace.Size(2) cellpathinfo.Dataspace.Size(3)]);

    for tp=firsttp:lasttp
        cellpath{tp} = cellpath_mat(:,:,tp);
    end
end

fid = H5F.open(fullfile(targetfolder,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,sisterList_name,'H5P_DEFAULT')
    H5F.close(fid);
    
    % Read in the sisterList matrix
    sisterList_mat = h5read(fullfile(targetfolder,H5filename),sisterList_name,[1 1 1], [sisterListinfo.Dataspace.Size(1) sisterListinfo.Dataspace.Size(2) sisterListinfo.Dataspace.Size(3)]);

    for tp=firsttp:lasttp
        sisterList{tp} = sisterList_mat(:,:,tp);
    end
end

fid = H5F.open(fullfile(targetfolder,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,bg_name,'H5P_DEFAULT')
    H5F.close(fid);
    
    % Read in the background point matrix, currently are not being used
    bg_mat = h5read(fullfile(targetfolder,H5filename),bg_name,[1 1 1], [bginfo.Dataspace.Size(1) bginfo.Dataspace.Size(2) bginfo.Dataspace.Size(3)]);

    for tp=firsttp:lasttp
        bg{tp} = bg_mat(:,:,tp);
    end
end


fid = H5F.open(fullfile(targetfolder,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');

if ~H5L.exists(fid,imagesegment,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' imagesegment]);
else
    H5L.delete(fid,imagesegment,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' imagesegment]);
    H5F.close(fid);
end

% Initialize image segment matrix
h5create(fullfile(targetfolder,H5filename), imagesegment, [size(cellpath{lasttp},1), lasttp, 3, 2*cellsize+1, 2*cellsize+1], 'Datatype', 'uint8', 'ChunkSize', [1, 1, 3, 2*cellsize+1, 2*cellsize+1], 'Deflate', 9);

selected_cells = [];
templateinfo = loadimageinfo(filetype,fileformat,[row col field plane templateCH],firsttp,channelnames,targetfolder);
imwidth = templateinfo.Width;
imheight = templateinfo.Height;

CellInd = find(cellpath{lasttp}(:,1)~=-1 & cellpath{lasttp}(:,2)~=-1)';
nucmask = cell(size(cellpath{lasttp},1),1);
cytomask = cell(size(cellpath{lasttp},1),1);
cellmask = cell(size(cellpath{lasttp},1),1);

% First convert cellpath coordinates to be independent of sisterList,
% meaning no cells will have coordinate of -1,-1 due to sistering
[new_cellpath,new_sisterList] = removeSister(cellpath,sisterList,firsttp,lasttp,1:length(cellpath{lasttp}));
tic
for cellNo = CellInd
    
    % Check to make sure that the current path is not out-of-bound in any
    % frames and contain the necessary coordinates
    [pathLogic,c_cellpath] = check_path(new_cellpath,cellNo,firsttp,lasttp,imheight,imwidth);
    
    if pathLogic
        % If good track, adding current cell to the selected_cells List
        selected_cells = [selected_cells;cellNo];
        
        template = zeros(2*cellsize+1,2*cellsize+1,lasttp);
        cellIm = zeros(2*cellsize+1,2*cellsize+1,lasttp);
        newsegments = zeros(firsttp,lasttp,3,size(template,1),size(template,2), 'uint8');
        display(['r' num2str(row) 'c' num2str(col) 'f' num2str(field) ',Processing cell:' num2str(cellNo)]);
        for tp = firsttp:lasttp
            if c_cellpath(tp,1)>0 && c_cellpath(tp,2)>0
                
                % Determine the image coordinate
                xL=max(c_cellpath(tp,1)-cellsize,1);
                xR=min(c_cellpath(tp,1)+cellsize,imwidth);
                yL=max(c_cellpath(tp,2)-cellsize,1);
                yR=min(c_cellpath(tp,2)+cellsize,imheight);
                
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
                
                neighborDistX = (new_cellpath{tp}(:,1)-c_cellpath(tp,1)); 
                neighborDistY = (new_cellpath{tp}(:,2)-c_cellpath(tp,2));
                neighborDist = (neighborDistX.^2+neighborDistY.^2).^(0.5);
                closebyCells = find(  neighborDist>0 & (neighborDist < cellsize*1.414) );
                cellCoord = [c_cellpath(tp,1)-xL c_cellpath(tp,2)-yL];
                neighborCoord = [new_cellpath{tp}(closebyCells,1)-xL new_cellpath{tp}(closebyCells,2)-yL];
                
                % Load the original images for the specified crop area, for
                % both template channel and whole cell channel
                template(borderY,borderX,tp) = loadimage(filetype,fileformat,[row col field plane templateCH -1],tp,channelnames,{[yL yR], [xL xR]},targetfolder);
                cellIm(borderY,borderX,tp) =   loadimage(filetype,fileformat,[row col field plane cellCH1 -1],tp,channelnames,{[yL yR], [xL xR]},targetfolder);
                
                % Determine nuclear, cytosolic and cell segment
                [nucMask,cellMask] = templateToCentroid(template(:,:,tp),cellIm(:,:,tp),cellCoord,neighborCoord,dilateSize,borderX,borderY);
                
                if ~isempty(find(nucMask==1,1))
                        newsegments(1,tp,1,:,:) = nucMask; %nuclei
                        newsegments(1,tp,2,:,:) = cellMask; %cell
                        newsegments(1,tp,3,:,:) = bwmorph(squeeze(newsegments(1,tp,2,:,:)),'erode',1) & ~squeeze(newsegments(1,tp,1,:,:)); % cytosol
                end
            end
            
        end
        t=toc;
        display(['r' num2str(row) 'c' num2str(col) 'f' num2str(field) ',Done Processing cell:' num2str(cellNo), ' Elapsed time:' num2str(t) 's']);
        % write the output segment to HDF5 file
        % order of the imagesegment matrix is : cell,time,segments,x,y
        h5write(fullfile(targetfolder,H5filename), imagesegment, newsegments, [cellNo firsttp 1 1 1], [1 lasttp-firsttp+1 3 size(template,1) size(template,2)]);

        %[uv sv] = memory;
        %fprintf('%f GB\n', uv.MemUsedMATLAB/1e9);
    end
end

% write the selectedcells list to HDF5 file
fid = H5F.open(fullfile(targetfolder,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');

if ~H5L.exists(fid,selectedcells_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' selectedcells_name]);
else
    H5L.delete(fid,selectedcells_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' selectedcells_name]);
    H5F.close(fid);
end

h5create(fullfile(targetfolder,H5filename), selectedcells_name, length(selected_cells), 'Datatype', 'uint32');
h5write(fullfile(targetfolder,H5filename), selectedcells_name, uint32(selected_cells));


% This code is used to convert cellTrack coordinates to individual separate
% coordinate without the embeded sister information
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

display('Done removing sisters');

function [pathLogic,ind_cellpath] = check_path(cellpath,cellNo,firsttp,lasttp,imheight,imwidth)
% check if this path is good
%#1 if acquiring and being inferior sister, combine trace of self with
%prior-sister
% has sister?

ind_cellpath = zeros(lasttp,2);

for t=firsttp:lasttp
    ind_cellpath(t,:) = cellpath{t}(cellNo,:);
end
deathInd = find(ind_cellpath(:,1)==-2,1,'first');

if isempty(deathInd)
    temp_cellpath = ind_cellpath;
else
    temp_cellpath = ind_cellpath(1:(deathInd-1),:);
end

% Containing large x-y drift in comparison to background points ?

% Is there any point that lies outside the image area.
if isempty(  find(temp_cellpath(:,1)>imwidth | temp_cellpath(:,1)<0 | temp_cellpath(:,2)>imheight | temp_cellpath(:,2)<0,1)  )

    diff_X = diff(temp_cellpath(:,1));
    diff_Y = diff(temp_cellpath(:,2));
    med_XY = median([diff_X;diff_Y]);
    std_XY = std([diff_X;diff_Y]);
    
    if isempty(  find(diff_X > (med_XY+30*std_XY) | diff_Y > (med_XY+30*std_XY) ,1)  )
        pathLogic = 1;
    else
        pathLogic = 0;
    end
    
else
    display(['cell ' num2str(cellNo) ':cell track out of bound']);
    pathLogic = 0;
end

% This code is used to code original images from the different image input type 
function outputim = loadimage(filetype,fileformat,imlocation,tp,channelnames,pixelsreg,targetfolder)
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
            outputim = imread(fullfile(targetfolder,filename),'PixelRegion',pixelsreg); 
        case 2
            outputim = imread(fullfile(targetfolder,fileformat),'Index',totalCH*(tp-1)+channel1,'PixelRegion',pixelsreg);
        case 3
            filename = sprintf(fileformat,channelnames{channel1},tp);
            outputim = imread(fullfile(targetfolder,filename),'PixelRegion',pixelsreg);
    end

    outputim = mat2gray(outputim);

else
    switch filetype
        case 1
            filename = sprintf(fileformat,row,col,field,plane,channel1,tp);
            nominim = imread(fullfile(targetfolder,filename),'PixelRegion',pixelsreg);
            filename = sprintf(fileformat,row,col,field,plane,channel2,tp);
            denomin = imread(fullfile(targetfolder,filename),'PixelRegion',pixelsreg);
        case 2
            nominim = imread(fullfile(targetfolder,fileformat),'Index',totalCH*(tp-1)+channel1,'PixelRegion',pixelsreg);
            denomin = imread(fullfile(targetfolder,fileformat),'Index',totalCH*(tp-1)+channel2,'PixelRegion',pixelsreg);
        case 3
            filename = sprintf(fileformat,channelnames{channel1},tp);
            nominim = imread(fullfile(targetfolder,filename),'PixelRegion',pixelsreg);
            filename = sprintf(fileformat,channelnames{channel2},tp);
            denomin = imread(fullfile(targetfolder,filename),'PixelRegion',pixelsreg);
    end
    outputim = mat2gray((im2double(nominim))./(im2double(denomin)));
end

function iminfo = loadimageinfo(filetype,fileformat,imlocation,tp,channelnames,targetfolder)
row = imlocation(1);
col = imlocation(2);
field = imlocation(3);
plane = imlocation(4);
channel = imlocation(5);
totalCH = length(channelnames);

switch filetype
    case 1
        filename = sprintf(fileformat,row,col,field,plane,channel,tp);
        iminfo = imfinfo(fullfile(targetfolder,filename));
    case 2
        iminfo = imfinfo(fullfile(targetfolder,fileformat),'Index',totalCH*(tp-1)+channel);
    case 3
        filename = sprintf(fileformat,channelnames{channel},tp);
        iminfo = imfinfo(fullfile(targetfolder,filename));
end


% for contour tracking analysis, currently not being used
function seg = modchenvese(I,mu,BW,maxIter,cellCoord,neighborCoord,nucMask)
I = imadjust(I);
mask = bwperim(BW);
imAdj = I-mask;
imOut=cat(3,max(imAdj,mask),imAdj,imAdj);
phi0 = bwdist(mask)-bwdist(1-mask)+im2double(mask)-.5;
force = eps;

%-- Display settings
figure(2);
subplot(2,2,1); imshow(imOut); title('Input Image');axis equal;
subplot(2,2,2); contour(flipud(phi0), [0 0], 'r','LineWidth',1); title('initial contour');axis equal;

%-- End Display original image and mask

%-- Main loop


n=1;
touching = 0;
while touching~=1 && n <= maxIter
    
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
    %indicator = checkstop(old,new,dt);
    
    NCoord =find(bwselect(BW,neighborCoord(:,1),neighborCoord(:,2)),1);
    % intermediate output
    if isempty(NCoord) 
        n=n+10;
        if (mod(n,500) == 0)
            subplot(2,2,3);
            showphi(I,phi0,n);
        end
        seg = (phi0<=0);
        seg = imfill(seg,'holes');
    else
        n=n+10;
        if (mod(n,500) == 0)
            subplot(2,2,3);
            showphi(I,phi0,n);
        end
        seg = (phi0<=0);
        seg = imfill(seg,'holes');
        NIm =bwselect(seg,neighborCoord(:,1),neighborCoord(:,2));
        CIm =bwselect(seg,cellCoord(1),cellCoord(2));
        if ~isempty(find(NIm & CIm,1));
            break;
        end
    end
    
end;
seg=bwselect(bwmorph(seg,'open'),cellCoord(1),cellCoord(2)) | nucMask;
subplot(2,2,4); imshow(seg);
% 
% %make mask from SDF
% seg = phi0<=0; %-- Get mask from levelset
% 
% midY = round(size(seg,1)/2);
% midX = round(size(seg,2)/2);
% 
% if seg(midY,midX) > 0
%     seg = bwselect(seg,midX,midY);
% else
%     seg = ~seg;
%     seg = bwselect(seg,midX,midY);
% end


%subplot(2,2,4); imshow(seg); title('Global Region-Based Segmentation');


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



function [notp stagePos stageName waveName] = readndfile(sourcefolder,filename)
% Search for number of string matches per line.
notp=-1;
stagePos = [];
stageName = [];
waveName = [];

fullfile(sourcefolder,filename)

if exist(fullfile(sourcefolder,filename),'file')
    fid = fopen(fullfile(sourcefolder,filename));
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
