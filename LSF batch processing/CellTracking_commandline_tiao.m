function CellTracking_commandline_tiao(filetype,SourceF,row,col,field,plane,nucCH,cellCH,tps,increment,fileformat,channelnames)
% Example usage: CellTracking_commandline(pwd,2,6,1,1,3,[1 16]);
% celltrackOUT    =  MAT file that contains initial track points
% tps             = [<first frame>   <last frame>]

load celltrackingparameters;

currentPath = pwd;
eval('cd ..');
addpath(genpath([pwd filesep 'ThirdParty']),'-end');
cd(currentPath);

%Calculate derived parameters
minNucArea = round(0.25*avgNucDiameter^2);
minNewNucArea = round(0.5*avgNucDiameter^2);
maxMAL = 5*avgNucDiameter;

%Define derection of tracking
firsttp = tps(1);
lasttp = tps(end);
switch increment
    case 1
        firstFrame = tps(1);
        endFrame   = tps(end)-1; 
    case -1
        firstFrame = tps(end);
        endFrame   = tps(1)+1; 
end


%Define H5File
H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
cellpath_name = ['/field' num2str(field) '/cellpath'];
sisterList_name = ['/field' num2str(field) '/sisterList'];
bg_name = ['/field' num2str(field) '/bg'];

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
    
    for tp=firsttp:cellpathinfo.Dataspace.Size(3)
        cellpath{tp} = cellpath_mat(:,:,tp);
    end
else
    cellpath = [];
end

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,sisterList_name,'H5P_DEFAULT')
    H5F.close(fid);
    sisterListinfo = h5info(fullfile(SourceF,H5filename), sisterList_name);
    sisterList_mat = h5read(fullfile(SourceF,H5filename),sisterList_name,[1 1 1], [sisterListinfo.Dataspace.Size(1) sisterListinfo.Dataspace.Size(2) sisterListinfo.Dataspace.Size(3)]);
    
    for tp=firsttp:sisterListinfo.Dataspace.Size(3)
        sisterList{tp} = sisterList_mat(:,:,tp);
    end
else
    sisterList = [];
end


if isempty(cellpath) || isempty(sisterList)  
    display([H5filename ' does not contain necessary variables.']);
    return
end

currCoords = cellpath{firstFrame};
currSisterList = sisterList{firstFrame};


% initiate lineage
if exist('lineageList','var')
    currLineage = lineageList{firstFrame};
else
    nCells = size(currCoords,1);
    currLineage = - ones(nCells,3);
    lineageList{firstFrame} = currLineage;
end


%% segment 1st frame and find nuclear centroids
% load in 1st frame
nucframe  = loadimage(filetype,fileformat,[row col field plane nucCH],firstFrame,channelnames,SourceF);
cellframe = loadimage(filetype,fileformat,[row col field plane cellCH],firstFrame,channelnames,SourceF);

% determine which ones are live tracks
nTracks=size(currCoords,1);
liveTrack=true(nTracks,1);
for c=1:nTracks
    % check whether it's a terminated track
    tmp=currLineage(c,:);
    tmp=sum(tmp==-2);
    if tmp == 3
        liveTrack(c)=false;
    end
end

% call function to segment nuclei
[bwNuc, ~] = detectSeededNuclei(nucframe, currCoords(liveTrack,:), avgNucDiameter);
% find seed point with no segmented nuclei and place a sphere there
nucCent = zeros(size(nucframe));
for i = 1 : size(currCoords,1)
    if liveTrack(i)
        nucCent(currCoords(i,2),currCoords(i,1)) = 1;
    end
end
nucCent = nucCent & ~bwNuc;
if max(nucCent(:)) > 0
    bwNucTmp = imdilate(nucCent, strel('disk',round(avgNucDiameter/2)));
    tmp = bwdist(bwNuc);
    tmp = tmp < 2;
    bwNucTmp = ~tmp & bwNucTmp;
    bwNuc = bwNuc | bwNucTmp;
end
% update the nuclear centroids to intensity weighted centroids
lNuc = bwlabel(bwNuc);
for c=1:size(currCoords,1)
    idx=lNuc(currCoords(c,2),currCoords(c,1));
    im=im2double(nucframe).*double(lNuc==idx);
    th=median(im(im>0));
    bw=im>=th;
    props=regionprops(bw,'Area');
    maxArea=max(cat(1,props.Area));
    bw=bwareaopen(bw,maxArea);
    props=regionprops(bw,nucframe,'WeightedCentroid');
    currCoords(c,:)=round(props(1).WeightedCentroid);
end
% write nuclear mask to file
bwCell=cellSegmentation(cellframe,bwNuc,thresParam,minAreaRatio,minCytosolWidth);

% output folder to store overlay image
overlayF = [SourceF filesep 'overlay'];
if exist(overlayF,'dir') == 0
    mkdir(overlayF);
end
% output folder to store nuclear masks
nucMaskF = [SourceF filesep 'nuclearMask'];
if exist(nucMaskF,'dir') == 0
    mkdir(nucMaskF);
end
% output folder to store nuclear masks
cellMaskF = [SourceF filesep 'cellMask'];
if exist(cellMaskF,'dir') == 0
    mkdir(cellMaskF);
end



% generate output image
nucCent = zeros(size(nucframe));
for i = 1 : size(currCoords,1)
    if liveTrack(i)
        nucCent(currCoords(i,2),currCoords(i,1)) = 1;
    end
end
p=imdilate(nucCent,strel('disk',2));
p3=bwperim(bwNuc);
p3=p3|bwperim(bwNuc&~p3);
p4=bwperim(bwCell);
p4=p4|bwperim(bwCell&~p4);
im=imadjust(im2double(nucframe));
im2=imadjust(im2double(cellframe));
imOut=im2uint8(cat(3,max(max(p3,p),p4),max(max(im2-p,p3),p4),max(im-p-p3,p4)));

% save output files
writeimage(imOut,filetype,fileformat,[row col field plane 1],firstFrame,channelnames,overlayF);
writeimage(bwNuc,filetype,fileformat,[row col field plane nucCH],firstFrame,channelnames,nucMaskF);
writeimage(bwCell,filetype,fileformat,[row col field plane cellCH],firstFrame,channelnames,cellMaskF);

%% loop over the remaining frames to segment and track nuclei
% compute parameter for the template matching function
load fftexecutiontimes
% calculate optimal conditions for fft analysis
fftw('planner', 'hybrid');
opt = detbestlength2(FFTrv,FFTiv,IFFTiv,2*[outersize outersize],2*[cellsize cellsize],1,1);
[ySize, xSize] = size(nucframe);
optWholeIm = detbestlength2(FFTrv,FFTiv,IFFTiv,2*[ySize ySize],2*[ySize-2*maxWholeImShift ySize-2*maxWholeImShift],1,1);

bwNuc_overlaps = zeros(size(nucframe));
bgOffset = zeros(max(firstFrame,endFrame+increment),2);
for t=firstFrame:increment:endFrame
    display([H5filename ':' 'T' num2str(t)]);
    % update previous frame
    prevCoords = currCoords;
    prevCoords2 = prevCoords;
    prevSisterList = currSisterList;
    prevLineage = currLineage;
    previousframe = nucframe;
    
    % load in current nuclear frame and cell frame
    tp=t+increment;
    nucframe =  loadimage(filetype,fileformat,[row col field plane nucCH],tp,channelnames,SourceF);
    cellframe = loadimage(filetype,fileformat,[row col field plane cellCH],tp,channelnames,SourceF);
    
    % calculate the whole image offset
    template=previousframe(maxWholeImShift+1:ySize-maxWholeImShift,maxWholeImShift+1:ySize-maxWholeImShift);
    testframe=nucframe(1:ySize,1:ySize);
    [x1,y1,maxVal] = corrMatching2(testframe, template, optWholeIm);
    xOffsetWholeIm = x1 - ySize/2;
    yOffsetWholeIm = y1 - ySize/2;
    prevCoords2(:,1) = prevCoords2(:,1) + xOffsetWholeIm;
    prevCoords2(:,2) = prevCoords2(:,2) + yOffsetWholeIm;
    
    bgOffset(tp,:) = [xOffsetWholeIm yOffsetWholeIm];
    
    
    % pad images with random background comparible to image background
    [ySize, xSize] = size(previousframe);
    [iMax,iWidth]=contrastValues(previousframe);
    previousPad=iMax+(rand(ySize+2*outersize+2*maxWholeImShift,xSize+2*outersize+2*maxWholeImShift)-0.5)*2*iWidth;
    previousPad(outersize+maxWholeImShift+1:ySize+outersize+maxWholeImShift,outersize+maxWholeImShift+1:xSize+outersize+maxWholeImShift)=double(previousframe);
    [iMax,iWidth]=contrastValues(nucframe);
    currentPad=iMax+(rand(ySize+2*outersize+2*maxWholeImShift,xSize+2*outersize+2*maxWholeImShift)-0.5)*2*iWidth;
    currentPad(outersize+maxWholeImShift+1:ySize+outersize+maxWholeImShift,outersize+maxWholeImShift+1:xSize+outersize+maxWholeImShift)=double(nucframe);
    prevCoordsPad = round(prevCoords) + outersize + maxWholeImShift;
    prevCoordsPad2 = round(prevCoords2) + outersize + maxWholeImShift;
    
    % call template matching function to update centroid locations

    nTracks=size(prevCoords,1);
    currCoords=zeros(size(prevCoords));
    liveTrack=true(nTracks,1);
    currSisterList = prevSisterList;
    
    for c=1:nTracks
        
        % check whether it's a terminated track
        tmp=prevLineage(c,:);
        tmp=sum(tmp==-2);
        if tmp == 3
            currCoords(c,:) = [-2 -2]; %prevCoords(c,:);
            liveTrack(c)=false;
        else
            xL=prevCoordsPad(c,1)-cellsize;
            xR=prevCoordsPad(c,1)+cellsize-1;
            yL=prevCoordsPad(c,2)-cellsize;
            yR=prevCoordsPad(c,2)+cellsize-1;
            template = previousPad(yL:yR,xL:xR);
            xL=prevCoordsPad2(c,1)-outersize;
            xR=prevCoordsPad2(c,1)+outersize-1;
            yL=prevCoordsPad2(c,2)-outersize;
            yR=prevCoordsPad2(c,2)+outersize-1;
            testframe = currentPad(yL:yR,xL:xR);
            [x1,y1,maxVal] = corrMatching2(testframe, template, opt);
            if isempty(x1) || maxVal<similarityThres
                x1 = prevCoordsPad2(c,1)-xL;
                y1 = prevCoordsPad2(c,2)-yL;
            end
            currCoords(c,:)=round([xL+x1 yL+y1]);
            currCoords(c,:) = currCoords(c,:) - outersize - maxWholeImShift;
            
            tmp=currCoords(c,1);
            if tmp > xSize || tmp < 1
                currCoords(c,:) = [-1 -1]; %prevCoords(c,:);
                liveTrack(c)=false;
            end
            
            tmp=currCoords(c,2);
            if tmp > ySize || tmp < 1
                currCoords(c,:) = [-1 -1]; %prevCoords(c,:);
                liveTrack(c)=false;
            end
            
        end
    end


    
    % call nuclear segmentation function
    [bwNuc, bwOther] = detectSeededNuclei(nucframe, currCoords(liveTrack,:), avgNucDiameter);
    % find seed point with no segmented nuclei
    lNucCent = zeros(size(nucframe));
    for i = 1 : size(currCoords,1)
        if liveTrack(i)
            lNucCent(currCoords(i,2),currCoords(i,1)) = i;
        end
    end
    nucCent = lNucCent > 0;
    % make sure the remaining nuclear centroids don't touch the existing
    % nuclei detected
    nucCent = nucCent & ~bwNuc;
    tmp=bwdist(nucCent);
    tmp=tmp<2;
    bwNuc = ~tmp & bwNuc;
    % for seed point with no segmented nuclei and search for unclaimed
    % object nearby
    if max(nucCent(:)) > 0 && max(bwOther(:)) > 0
        distMat = bwdist(bwOther);
        idx = find(nucCent);
        lList = lNucCent(idx);
        distVal = distMat(idx);
        [minDistVal, minDistIdx] = min(distVal);
        while minDistVal <= maxNucMaskShift && max(nucCent(:)) > 0
            idx = idx(minDistIdx);
            lList = lList(minDistIdx);
            lOther = bwlabel(bwOther);
            nucCentTmp = zeros(size(nucframe));
            nucCentTmp(idx) = 1;
            tmp = bwdist(nucCentTmp);
            tmp = tmp <= minDistVal;
            tmp = unique(lOther(tmp));
            tmp(tmp == 0) = [];
            tmp = tmp(1);
            props = regionprops(lOther == tmp, nucframe, 'WeightedCentroid');
            currCoords(lList,:)=round(props(1).WeightedCentroid);
            bwNuc = bwNuc | lOther == tmp;
            bwOther = bwOther & ~bwNuc;
            nucCent(idx) = false;
            distMat = bwdist(bwOther);
            idx = find(nucCent);
            lList = lNucCent(idx);
            distVal = distMat(idx);
            [minDistVal, minDistIdx] = min(distVal);
            if isempty(minDistVal)
                minDistVal = maxNucMaskShift + 1;
            end
        end
    end
    % place a sphere at the centroids with no matching segmented nuclei
    bwNucMissing = false(size(bwNuc));
    if max(nucCent(:)) > 0
        tmp=bwdist(nucCent);
        tmp=imimposemin(tmp,nucCent);
        l=watershed(tmp);
        bwNucMissing = imdilate(nucCent, strel('disk',round(avgNucDiameter/2)));
        tmp = bwdist(bwNuc);
        tmp = tmp < 2;
        bwNucMissing = (~tmp & bwNucMissing) & l > 0;
        l = bwlabel(bwNucMissing);
        idx = unique(l(nucCent));
        idx(idx == 0) = [];
        bwNucMissing = ismember(l,idx);
        bwNuc = bwNuc | bwNucMissing;
        bwOther = bwOther & ~bwNuc;
    end
    % update nuclear centroids to intensity weighted centroids
    lNuc = bwlabel(bwNuc);
    nCent = size(currCoords,1);
    currLineage = zeros(nCent,3);
    currLineage(:,2:3) = prevLineage(:,1:2);
    for c = 1 : nCent
        % update current lineage list
        if ~liveTrack(c)
            currLineage(c,1) = -2;  % terminated tracks
        else
            if bwNucMissing(currCoords(c,2),currCoords(c,1))
                currLineage(c,1) = -2;  % missing nuclei is label as -2
            else
                currLineage(c,1) = c;
            end
            idx=lNuc(currCoords(c,2),currCoords(c,1));
            im=im2double(nucframe).*double(lNuc==idx);
            th=median(im(im>0));
            bw=im>=th;
            props=regionprops(bw,'Area');
            maxArea=max(cat(1,props.Area));
            bw=bwareaopen(bw,maxArea);
            props=regionprops(bw,nucframe,'WeightedCentroid');
            currCoords(c,:)=round(props(1).WeightedCentroid);
        end
    end

% --  identify new nuclei
%     bwOtherTmp = blobSegmentThresholdGeneral(nucframe, 'rosin', 1, 1, minNucArea);
%     lOther = bwlabel(bwOther);
%     idx = unique(lOther(bwOtherTmp));
%     idx(idx == 0) = [];
%     bwNucTmp = ismember(lOther,idx);
    bwNucTmp = bwOther;
    if max(bwNucTmp(:)) > 0
        [lNucTmp, nNucTmp] = bwlabel(bwNucTmp);
        S = regionprops(lNucTmp,'Perimeter','Area','MajorAxisLength');
        for n = 1: nNucTmp
            if S(n).MajorAxisLength > maxMAL || S(n).Area < minNewNucArea || 4*pi*S(n).Area/S(n).Perimeter^2 < ffactorCutoff
                bwNucTmp(lNucTmp == n) = false;
            end
        end
    end
    if max(bwNucTmp(:)) > 0
        [lNucTmp, nNucTmp] = bwlabel(bwNucTmp);
        for c = nCent+1 : nCent + nNucTmp
            im=im2double(nucframe).*double(lNucTmp == c-nCent);
            th=median(im(im>0));
            bw=im>=th;
            props=regionprops(bw,'Area');
            maxArea=max(cat(1,props.Area));
            bw=bwareaopen(bw,maxArea);
            props=regionprops(bw,nucframe,'WeightedCentroid');
            currCoords(c,:)=round(props(1).WeightedCentroid);
            % determine whether it's a cell moving into the field of view
            % or newly divided one
            if currCoords(c,1) < 1.5*avgNucDiameter || currCoords(c,1) > (xSize - 1.5*avgNucDiameter) || ...
                    currCoords(c,2) < 1.5*avgNucDiameter || currCoords(c,2) > (ySize - 1.5*avgNucDiameter)
                currLineage(c,1:3)=-1;  % most likely moved into the field of view
                currSisterList(c,1:end) = -1;
            else
                % most likely a divided one, so find the closest nucleus as
                % its sister
                tmp = (currCoords(1:nCent,1) - currCoords(c,1)).^2 + (currCoords(1:nCent,2) - currCoords(c,2)).^2;
                tmp(~liveTrack)=max(tmp);
                [~,idx] = min(tmp);
                currLineage(c,1:3) = currLineage(idx,1:3);
                % Update siterList of parent cells
                nextSisInd = find(currSisterList(idx,:)==-1,1,'first');
                currSisterList(idx,nextSisInd) = c;
                
                % Update sisterList of newly created cells
                currSisterList(c,1) = idx;
                currSisterList(c,2:end) = -1;
                
            end
        end
        bwNuc = bwNuc | bwNucTmp;
    end
    
    % update cellSegment
    bwCell=cellSegmentation(cellframe,bwNuc,thresParam,minAreaRatio,minCytosolWidth);
    
    cellpath{tp} = currCoords;
    sisterList{tp} = currSisterList;
    lineageList{tp} = currLineage;
    
    if increment > 0 && length(cellpath{tp-increment}) < length(cellpath{tp})
        for t2=firstFrame:increment:tp-increment
            for c = length(cellpath{tp-increment})+1 : length(cellpath{tp})
                cellpath{t2}(c,:) = [-1 -1];
                sisterList{t2}(c,:) = -1 * ones(1,size(sisterList{tp},2));
            end
        end
    end

    %create output image
    nucCent=zeros(size(nucframe));
    for i=1:nCent
        if liveTrack(i)
            nucCent(currCoords(i,2),currCoords(i,1))=1;
        end
    end
    p=imdilate(nucCent,strel('disk',2));
    nucCent=zeros(size(nucframe));
    for i=nCent+1:size(currCoords,1)
        nucCent(currCoords(i,2),currCoords(i,1))=1;
    end
    
    bwNuc_overlaps = bwNuc_overlaps | bwNuc;
    
    %For showing result
    p2=imdilate(nucCent,strel('disk',2));
    p3=bwperim(bwNuc);
    p3=p3|bwperim(bwNuc&~p3);
    p4=bwperim(bwCell);
    p4=p4|bwperim(bwCell&~p4);
    im=imadjust(im2double(nucframe));
    im2=imadjust(im2double(cellframe));
    imOut=im2uint8(cat(3,max(max(p3-p2,p),p4),max(max(im2-p-p2,p3),p4),max(max(im-p-p3,p2),p4)));
    % write output images
    writeimage(imOut,filetype,fileformat,[row col field plane 1],tp,channelnames,overlayF);
    writeimage(bwNuc,filetype,fileformat,[row col field plane nucCH],tp,channelnames,nucMaskF);
    writeimage(bwCell,filetype,fileformat,[row col field plane cellCH],tp,channelnames,cellMaskF);
    clear imOut bwNuc bwCell p2 p3 p4 im im2
end

bg_im = ~imfill(imdilate(bwNuc_overlaps,strel('disk',60)),'holes');
bg_im(1:end,[1:30 size(bg_im,2)-30:size(bg_im,2)]) = 0;
bg_im([1:30 size(bg_im,1)-30:size(bg_im,1)],1:end) = 0;
bg_im = imopen(bg_im,strel('disk',40));
bg_im_p = bwmorph(bg_im,'shrink',Inf);

[bg_Y,bg_X] = find(bg_im_p);

% Saving results
cellpath_mat = -1*(ones(size(cellpath{tp},1),2,length(cellpath)));
sisterList_mat = -1*(ones(size(sisterList{tp},1),size(sisterList{tp},2),length(sisterList)));
bg_mat = [];
for s_tp=firstFrame:increment:(endFrame+increment)
    if ~isempty(cellpath{s_tp})
        cellpath_mat(:,:,s_tp) = cellpath{s_tp};
        sisterList_mat(:,:,s_tp) = sisterList{s_tp};
        offset = bgOffset(s_tp,:);
        if ~isempty(bg_X) || ~isempty(bg_Y)
            
            bg_mat(:,:,s_tp) = round([bg_X bg_Y] +[offset(1)*ones(size(bg_X)) offset(2)*ones(size(bg_Y))]);
        else
            
            bg_mat(:,:,s_tp) = round([offset(1) offset(2)]);
        end
    end
end

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' cellpath_name]);
else
    H5L.delete(fid,cellpath_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' cellpath_name]);
    H5F.close(fid);
end
h5create(fullfile(SourceF,H5filename), cellpath_name, [size(cellpath_mat,1), size(cellpath_mat,2), size(cellpath_mat,3)], 'Datatype', 'double', 'ChunkSize', [1, size(cellpath_mat,2), size(cellpath_mat,3)], 'Deflate', 9);
h5write(fullfile(SourceF,H5filename), cellpath_name, cellpath_mat, [1 1 1], [size(cellpath_mat,1) size(cellpath_mat,2) size(cellpath_mat,3)]);

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,sisterList_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' sisterList_name]);
else
    H5L.delete(fid,sisterList_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' sisterList_name]);
    H5F.close(fid);
end
h5create(fullfile(SourceF,H5filename), sisterList_name, [size(sisterList_mat,1), size(sisterList_mat,2), size(sisterList_mat,3)], 'Datatype', 'double', 'ChunkSize', [1, size(sisterList_mat,2), size(sisterList_mat,3)], 'Deflate', 9);
h5write(fullfile(SourceF,H5filename), sisterList_name, sisterList_mat, [1 1 1], [size(sisterList_mat,1) size(sisterList_mat,2) size(sisterList_mat,3)]);

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,bg_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' bg_name]);
else
    H5L.delete(fid,bg_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' bg_name]);
    H5F.close(fid);
end
h5create(fullfile(SourceF,H5filename), bg_name, [size(bg_mat,1), size(bg_mat,2), size(bg_mat,3)], 'Datatype', 'double', 'ChunkSize', [1, size(bg_mat,2), size(bg_mat,3)], 'Deflate', 9);
h5write(fullfile(SourceF,H5filename), bg_name, bg_mat, [1 1 1], [size(bg_mat,1) size(bg_mat,2) size(bg_mat,3)]);


% Finished all processes
display(['Successfully tracked frame ' num2str(tps(1)) ':' num2str(tps(end)) ' of ' H5filename]);


function bwCell=cellSegmentation(cellIm,bwNuc,thresParam,minAreaRatio,minCytosolWidth)

    
%% cell segmentation
% bwNuc=bwNuc|bwOther;
% cellIm=imread(['02032013-r1_w1Camera-YFP_s1_t' num2str(n) '.TIF']);
cellIm=im2double(cellIm);
cellImFilled=cellIm;
cellImFilled(bwNuc)=1;
cellImS=filterGauss2D(cellIm,1);

% dilate nuclei by the minCytosolWidth pixels
l=-bwdist(~bwNuc);
l=imimposemin(l,bwNuc);
l=watershed(l);
nucDiv=l==0;
dilateNuc=bwdist(bwNuc);
dilateNuc=dilateNuc<=minCytosolWidth;
dilateNuc=dilateNuc&~nucDiv;
% used dilated nuclei to plug to holes in the whole cells channle, then
% perform watershed segmentation
tmp=-filterGauss2D(cellImFilled,1);
%     tmp(dilateNuc)=-1;
l=tmp;
l=imimposemin(l,bwNuc);
l=watershed(l);
% threshold whole cell nucCH
[iMax,iWidth]=contrastValues(im2uint16(cellImS));
thCell=(iMax+thresParam*iWidth)/65535;
bwCell=cellImS>thCell;
bwCell=bwCell|dilateNuc;
bwCell=bwCell&l>0;
% remove pathes w/ no corresponding nuclei
lCell=bwlabel(bwCell);
idx=unique(lCell(bwNuc));
idx(idx==0)=[];
bwCell=ismember(lCell,idx);

    % refine threshold for mitotic cells
    [lCell,nCell]=bwlabel(bwCell);
    lCell=double(lCell);
    lNuc=lCell.*double(bwNuc);
    th=ones(size(cellIm));
    for i=1:nCell
        tmp=cellImS(lNuc==i);
        medInt=median(tmp);
        stdInt=std(double(tmp));
        th(lCell==i)=medInt-stdInt;
    end
    bw=cellImS>th;
    bw=bw|bwNuc;
    l=lCell.*double(bw);
    props=regionprops(lNuc,'Area');
    nucArea=cat(1,props.Area);
    props=regionprops(l,'Area');
    bwArea=cat(1,props.Area);
    areaRatio=bwArea./nucArea;
    idx=find(areaRatio < minAreaRatio);

    bwMitoCell=false(size(cellIm));
    for i=1:numel(idx)
        tmp=lNuc==idx(i);
        tmp=bwdist(tmp);
        tmp=tmp<=minCytosolWidth;
        tmp=tmp&(lCell==idx(i));
        bwMitoCell=bwMitoCell|tmp;
        bwCell(lCell == idx(i)) = false;
    end
    bwCell = bwCell | bwMitoCell;


function outputim = loadimage(filetype,fileformat,imlocation,tp,channelnames,SourceF)
row = imlocation(1);
col = imlocation(2);
field = imlocation(3);
plane = imlocation(4);
CH = imlocation(5);
totalCH = length(channelnames);
outputim = [];

switch filetype
    case 1
        
        filename = sprintf(fileformat,row,col,field,plane,CH,tp);
        if exist(fullfile(SourceF,filename),'file')
            outputim = imread(fullfile(SourceF,filename));
        end
    case 2
        if exist(fileformat,'file')
            outputim = imread(fileformat,'Index',totalCH*(tp-1)+CH);
        end
    case 3
        filename = sprintf(fileformat,channelnames{CH},tp);
        if exist(fullfile(SourceF,filename),'file');
            outputim = imread(fullfile(SourceF,filename));
        end
        
end

function writeimage(im,filetype,fileformat,imlocation,tp,channelnames,destF)
row = imlocation(1);
col = imlocation(2);
field = imlocation(3);
plane = imlocation(4);
CH = imlocation(5);

switch filetype
    case 1
        
        filename = sprintf(fileformat,row,col,field,plane,CH,tp);
        imwrite(im,fullfile(destF,filename),'tif','compression','none');
    case 2
        
    case 3
        filename = sprintf(fileformat,channelnames{CH},tp);
        imwrite(im,fullfile(destF,filename),'tif','compression','none');
        
end

function [x y BW] = templateToCentroid(M,xg,yg)
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

if isempty(find(BW==0)) | isempty(find(BW==1))
    x = xg;
    y = yg;

else
    x = round(S.Centroid(1));
    y = round(S.Centroid(2));
end


function [out]=fftolamopt2(a,b,opt,shape)
% [out]=fftolamopt2(a,b,siz1,siz2,shape)
%
% Overlap-add method FFT-based 2D convolution
% Example:
%   load fftexecutiontimes;                                                        % load FFTrv, FFTiv and IFFTiv in workspace
%   a   = rand(500,500);                                                           % first image
%   b   = rand(340,220);                                                           % second image
%   opt = detbestlength2(FFTrv,FFTiv,IFFTiv,size(a),size(b),isreal(a),isreal(b));  % optimized parameters
%   y0  = fftolamopt2(a,b,opt);                                                    % equivalent to y0 = conv2(a,b);
%
% INPUT
% a:     first image (2D double matrix)
% b:     second image (2D double matrix)
% opt:   the optimized parameters calculated by detbestlength.m function
%        opt = detbestlength(FFTrv,FFTiv,IFFTiv,size(a),size(b));
% shape: returns a subsection of the 2D convolution with size specified by
%        'shape':
%          'full'  - (default) returns the full 2-D convolution,
%          'same'  - returns the central part of the convolution
%                    that is the same size as A.
%          'valid' - returns only those parts of the convolution
%                    that are computed without the zero-padded
%                    edges. size(C) = [ma-mb+1,na-nb+1] when
%                    all(size(A) >= size(B)), otherwise C is empty.
% See also conv2.
% OUTPUT
% out:   2D convolution of a and b matrices: out = conv2(a,b);


% Original size
[z1x,z1y] = size(a);
[z2x,z2y] = size(b);

% Reverse a and b if necessary
if opt.inverse
    atemp = a;
    a     = b;
    b     = atemp;
end

fftorder  = zeros(2,1);
ifftorder = zeros(2,1);
fftsize   = zeros(2,1);
filterord = zeros(2,1);
filtersiz = zeros(2,1);

if (opt.fftxfirst == 1)
    fftorder(1)  = 1;
    fftorder(2)  = 2;
    fftsize(1)   = opt.nfftx;
    fftsize(2)   = opt.nffty;
else
    fftorder(1)  = 2;
    fftorder(2)  = 1;
    fftsize(1)   = opt.nffty;
    fftsize(2)   = opt.nfftx;
end


if (opt.ifftxfirst == 1)
    ifftorder(1) = 1;
    ifftorder(2) = 2;
else
    ifftorder(1) = 2;
    ifftorder(2) = 1;
end

if opt.filterxfirst==1
    filterord(1) = 1;
    filterord(2) = 2;

    filtersiz(1) = opt.nfftx;
    filtersiz(2) = opt.nffty;
else
    filterord(1) = 2;
    filterord(2) = 1;

    filtersiz(1) = opt.nffty;
    filtersiz(2) = opt.nfftx;
end

siz1          = opt.nfftx;
siz2          = opt.nffty;

[ax,ay]       = size(a);
[bx,by]       = size(b);
dimx          = ax+bx-1;
dimy          = ay+by-1;
nfftx         = siz1;
nffty         = siz2;
Lx            = nfftx-bx+1;
Ly            = nffty-by+1;
B             = fft(fft(b,filtersiz(1),filterord(1)),filtersiz(2),filterord(2));
out           = zeros(dimx,dimy);
x0 = 1;
while x0 <= ax
    x1   = min(x0+Lx-1,ax);
    y0   = 1;
    endx = min(dimx,x0+nfftx-1);
    while y0 <= ay
        y1                   = min(y0+Ly-1,ay);
        endy                 = min(dimy,y0+nffty-1);
        X                    = fft(fft(a(x0:x1,y0:y1),fftsize(1),fftorder(1)),fftsize(2),fftorder(2));
        Y                    = ifft(ifft(X.*B,[],ifftorder(1)),[],ifftorder(2));
        out(x0:endx,y0:endy) = out(x0:endx,y0:endy)+Y(1:(endx-x0+1),1:(endy-y0+1));
        y0                   = y0+Ly;
    end
    x0 = x0+Lx;
end
if isreal(a) && isreal(b)
    out=real(out);
end
if nargin<4 || strcmp(shape,'full')
    return;
end
if strcmp(shape,'valid')
    if ((z1x<z2x)||(z1y<z2y))
        out = [];
    else
        px  = z2x;
        py  = z2y;
        out = out(px:px+z1x-z2x,py:py+z1y-z2y);
    end
    return;
end
if strcmp(shape,'same')
    px  = ((z2x-1)+mod((z2x-1),2))/2;
    py  = ((z2y-1)+mod((z2y-1),2))/2;
    out = out(px+1:px+z1x,py+1:py+z1y);
    return;
end

function [out]=detbestlength2(FFTrv,FFTiv,IFFTiv,size1,size2,isreal1,isreal2)
% [out]=detbestlength2(FFTrv,FFTiv,IFFTiv,size1,size2,isreal1,isreal2)
% Determine the best parameters for Overlap-Add FFT-based convolution.
%
% INPUT
% FFTrv:   vector with costs of FFT for real 1d vectors
% FFTiv:   vector with costs of FFT for complex 1d vectors
% IFFTiv:  vector with costs of IFFT for complex 1d vectors
% size1:   size(first_image)
% size2:   size(second_image)
% isreal1: 1 if first image is real, 0 otherwise (complex)
% isreal2: 1 if second image is real, 0 otherwise (complex)
% OUTPUT
% out:    the optimized parameters:
%         out.inverse:     if 1 the two input have to be inverted
%         out.fftxfirst:   if one the image has to be fft first along
%                          x-dimension
%         out.ifftxfirst:  if one the product of spectra has to be ifft
%                          first along x-dimensio
%         out.nfftx:       the best length for fft transform along
%                          x-dimension
%         out.nffty:       the best length for fft transform along
%                          y-dimension
%         out.filterxfirst if 1 the filter has to be fft fisrt alng
%                          x-dimension
%

out           = [];
% the 3 input vectors have to be the same length
L             = length(FFTrv);
% a default value (just as Inf)
infinitevalue = 99*10^99;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----------------------------------------------------- a image and b filter
if isreal1 && isreal2
    ax = size1(1);
    ay = size1(2);
    bx = size2(1);
    by = size2(2);

    val0 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));

                    if cv1<val0
                        val0 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val0
                        val0 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val0
                        val0 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val0
                        val0 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val0
                        val0 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val0
                        val0 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val0
                        val0 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val0
                        val0 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %----------------------------------------------------- a filter and b image
    ax = size2(1);
    ay = size2(2);
    bx = size1(1);
    by = size1(2);


    val1 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));

                    if cv1<val1
                        val1 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val1
                        val1 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val1
                        val1 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val1
                        val1 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val1
                        val1 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val1
                        val1 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val1
                        val1 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val1
                        val1 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %--------------------------------------------------------------------------
    if val1<val0
        out.inverse = 1;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    else
        out.inverse = 0;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    end
    out.nfftx = x;
    out.nffty = y;
    if t==1
        out.filterxfirst = 0;
    else
        out.filterxfirst = 1;
    end
    return;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----------------------------------------------------- a image and b filter
if ~isreal1 && ~isreal2
    ax = size1(1);
    ay = size1(2);
    bx = size2(1);
    by = size2(2);

    val0 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));

                    if cv1<val0
                        val0 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val0
                        val0 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val0
                        val0 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val0
                        val0 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val0
                        val0 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val0
                        val0 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val0
                        val0 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val0
                        val0 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %----------------------------------------------------- a filter and b image
    ax = size2(1);
    ay = size2(2);
    bx = size1(1);
    by = size1(2);


    val1 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTiv(ii)~=0 && FFTiv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));

                    if cv1<val1
                        val1 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val1
                        val1 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val1
                        val1 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val1
                        val1 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val1
                        val1 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val1
                        val1 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val1
                        val1 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val1
                        val1 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %--------------------------------------------------------------------------
    if val1<val0
        out.inverse = 1;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    else
        out.inverse = 0;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    end
    out.nfftx = x;
    out.nffty = y;
    if t==1
        out.filterxfirst = 0;
    else
        out.filterxfirst = 1;
    end
    return;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----------------------------------------------------- a image and b filter
if isreal1 && ~isreal2
    ax = size1(1);
    ay = size1(2);
    bx = size2(1);
    by = size2(2);

    val0 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));

                    if cv1<val0
                        val0 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val0
                        val0 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val0
                        val0 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val0
                        val0 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val0
                        val0 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val0
                        val0 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val0
                        val0 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val0
                        val0 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %----------------------------------------------------- a filter and b image
    ax = size2(1);
    ay = size2(2);
    bx = size1(1);
    by = size1(2);


    val1 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));

                    if cv1<val1
                        val1 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val1
                        val1 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val1
                        val1 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val1
                        val1 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val1
                        val1 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val1
                        val1 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val1
                        val1 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val1
                        val1 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %--------------------------------------------------------------------------
    if val1<val0
        out.inverse = 1;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    else
        out.inverse = 0;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    end
    out.nfftx = x;
    out.nffty = y;
    if t==1
        out.filterxfirst = 0;
    else
        out.filterxfirst = 1;
    end
    return;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----------------------------------------------------- a image and b filter
if ~isreal1 && isreal2
    ax = size1(1);
    ay = size1(2);
    bx = size2(1);
    by = size2(2);

    val0 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));

                    if cv1<val0
                        val0 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val0
                        val0 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val0
                        val0 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val0
                        val0 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val0
                        val0 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val0
                        val0 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val0
                        val0 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val0
                        val0 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %----------------------------------------------------- a filter and b image
    ax = size2(1);
    ay = size2(2);
    bx = size1(1);
    by = size1(2);


    val1 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));

                    if cv1<val1
                        val1 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val1
                        val1 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val1
                        val1 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val1
                        val1 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val1
                        val1 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val1
                        val1 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val1
                        val1 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val1
                        val1 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %--------------------------------------------------------------------------
    if val1<val0
        out.inverse = 1;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    else
        out.inverse = 0;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    end
    out.nfftx = x;
    out.nffty = y;
    if t==1
        out.filterxfirst = 0;
    else
        out.filterxfirst = 1;
    end
    return;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



function [x y maxVal] = corrMatching2(frameImg, templateImg,opt)
% -------------------------------------------------------------------------
% Function corrMatching: Template Matching using Correlation Coefficients
% Inputs: 
%           frameImg = gray or color frame image
%           templateImg = gray or color template image
%           threshC = threshold of rejecting detected region (default = .75)
%                     e.g. if the detected region has a corrCoef>threshC
%                     then the algorithm accepts it as a detection,
%                     otherwise rejects it as a false alarm.
% Output: 
%           corrScore = 2D matrix of correlation coefficients
%           boundingBox = [upperLeftPixel.y upperLeftPixel.x height width]
%
% -------------------------------------------------------------------------
% By Yue Wu (Rex)
% Department of Electrical and Computer Engineering
% Tufts University
% Medford, MA
% 08/30/2010
% -------------------------------------------------------------------------

% 1. initialization

if size(frameImg,3) ~=1
    frameGray = rgb2gray(frameImg);
else
    frameGray = frameImg;
end
frameGray = double(frameGray);

if size(templateImg,3) ~=1
    templateGray = rgb2gray(templateImg);
else
    templateGray = templateImg;
end

templateGray = double(templateGray);
[templateHeight,templateWidth] = size(templateGray);


% 2. correlation calculation

frameMean = fftolamopt2(frameGray,ones(size(templateGray))./numel(templateGray),opt);
templateMean = mean(templateGray(:));

corrPartI = fftolamopt2(frameGray,fliplr(flipud(templateGray-templateMean)),opt)./numel(templateGray);
corrPartII = frameMean.*sum(templateGray(:)-templateMean);
stdFrame = sqrt(fftolamopt2(frameGray.^2,ones(size(templateGray))./numel(templateGray),opt)-frameMean.^2);
stdTemplate = std(templateGray(:));
corrScore = (corrPartI-corrPartII)./(stdFrame.*stdTemplate);
% figure(3);imshow(corrScore);
% 3. finding most likely region

[maxVal,maxIdx] = max(corrScore(:));
[y, x] = ind2sub([size(corrScore,1),size(corrScore,2)],maxIdx);
%figure(3); hold on; plot(x,y,'rx');
x = x-round(templateWidth/2);
y = y-round(templateHeight/2);

