function CellSeeding_commandline(filetype,SourceF,row,col,field,plane,channel,tps,increment,forceseeding,minNucDiameter,maxNucDiameter,fileformat,channelnames)


currentPath = pwd;
eval('cd ..');
addpath(genpath([pwd filesep 'ThirdParty']),'-end');
cd(currentPath);
switch increment
    case 1
        firstFrame = tps(1);
    case -1
        firstFrame = tps(end);
end

H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
cellpath_name = ['/field' num2str(field) '/cellpath'];

% output folder to store overlay image
overlayF = [SourceF filesep 'overlay'];
if exist(overlayF,'dir') == 0
    mkdir(overlayF);
end

if exist(fullfile(SourceF,H5filename),'file') && ~forceseeding
    display([H5filename 'existed and will not be modified.']);
    return
elseif exist(fullfile(SourceF,H5filename),'file') && forceseeding
    fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if ~H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
        H5F.close(fid);
        display(['Initializing ' H5filename ':' cellpath_name]);
    else
        H5L.delete(fid,cellpath_name,'H5P_DEFAULT');
        display(['Overwriting ' H5filename ':' cellpath_name]);
        H5F.close(fid);

    end
end

firstframe = loadimage(filetype,fileformat,[row col field plane channel],firstFrame,channelnames,SourceF);
padWidth =10;
minFormfactor=0.7;
[ySize, xSize] = size(firstframe);
firstframe=im2uint16(firstframe);
[iMax, iWidth] = contrastValues(firstframe);
tmp = iMax + (rand(ySize+2*padWidth, xSize+2*padWidth)-0.5) * 2*iWidth;
tmp(padWidth+1:ySize+padWidth, padWidth+1:xSize+padWidth) = firstframe;
firstframe = tmp/max(tmp(:));
maxNucArea=round(pi*maxNucDiameter^2/4);

bwNuc = false(size(tmp));
testSet = linspace(1.1,4,5);

for s = 1:length(testSet)
    
    [im2,~]  = edge(tmp,'canny',0,testSet(s));
    im3 = ~im2;
    bw=im3-bwareaopen(im3,maxNucArea,4);
    bw=imfill(bw,'holes');
    CC = bwconncomp(bw, 8);
    L = labelmatrix(CC);
    S = regionprops(CC,'EquivDiameter','Area','Perimeter');
    nucArea=cat(1,S.Area);
    nucPerim=cat(1,S.Perimeter);
    nucEquiDiameter=cat(1,S.EquivDiameter);
    nucFormfactor=4*pi*nucArea./(nucPerim.^2);
    
    GoodCells = find(nucEquiDiameter >= minNucDiameter & ...
        nucEquiDiameter < maxNucDiameter & ...
        nucFormfactor>minFormfactor);
    bwNuc = bwNuc | ismember(L,GoodCells);
end

% remove padded edge region
bwNuc = bwNuc(padWidth+1:ySize+padWidth, padWidth+1:xSize+padWidth);
firstframe = firstframe(padWidth+1:ySize+padWidth, padWidth+1:xSize+padWidth);

% determine centroid
CC = bwconncomp(bwNuc, 8);
S = regionprops(CC,'Centroid');
numCell = length(S);

if numCell>0
    cellpath_mat = -1*(ones(numCell,2,tps(end)));
    nucCent=zeros(size(firstframe));
    for c=1:numCell
        cellpath_mat(c,1,firstFrame) = round(S(c).Centroid(1));
        cellpath_mat(c,2,firstFrame) = round(S(c).Centroid(2));
        nucCent(cellpath_mat(c,2,firstFrame),cellpath_mat(c,1,firstFrame))=1;
    end
    
    h5create(fullfile(SourceF,H5filename), cellpath_name, ...
        [size(cellpath_mat,1), size(cellpath_mat,2), size(cellpath_mat,3)], ...
        'Datatype', 'double', 'ChunkSize', [1, size(cellpath_mat,2), size(cellpath_mat,3)], 'Deflate', 9);
    h5write(fullfile(SourceF,H5filename), cellpath_name, cellpath_mat, [1 1 1], ...
        [size(cellpath_mat,1) size(cellpath_mat,2) size(cellpath_mat,3)]);
    

    p1=imdilate(nucCent,strel('disk',2));
    p2 = bwperim(bwNuc);
    p2 = bwmorph(p2,'dilate',1);
    im=imadjust(im2double(firstframe));
    imOut=cat(3,max(im-p2-p1,p1),max(im-p2-p1,p2),im-p2-p1 );
    writeimage(imOut,filetype,fileformat,[row col field plane channel],firstFrame,channelnames,overlayF);
    display(['Successfully seeded frame ' num2str(firstFrame) ' and saved in ' H5filename]);
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
        if exist(fileformat,'file')
            outputim = imread(fileformat,'Index',totalCH*(tp-1)+channel);
        end
    case 3
        filename = sprintf(fileformat,channelnames{channel},tp);
        if exist(fullfile(SourceF,filename),'file');
            outputim = imread(fullfile(SourceF,filename));
        end
end


