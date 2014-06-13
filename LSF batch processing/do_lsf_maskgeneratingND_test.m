function do_lsf_maskgeneratingND_test()
% Define information about input images and necessary parameters-----------
templateCH = 1;
cellsize = 18;
celldilatesize = 4;
%---------------------------------------
ndfilename ='0223201410A1806-r1.nd';
sourcefolder = 'Q:\sorger\data\NIC\Pat\03-23-2014-10A1806';
%------------------------------------------------

prefix = ndfilename(1:(end-3));
[notp,stagePos,stageName,channelnames] = readndfile(sourcefolder,ndfilename);
tps = [1 notp];
sites = [45 64];

if matlabpool('size') == 0
  matlabpool open;
end

parfor i=1:length(sites)
    site = sites(i);
    
    fileformat = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
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
    plane = 1;
    
    maskgen_commandline(3,sourcefolder,row,col,field,plane,templateCH,tps,fileformat,channelnames,cellsize,celldilatesize);
end



function maskgen_commandline(filetype,SourceF,row,col,field,plane,templateCH,tps,fileformat,channelnames,cellsize,celldilatesize)
firsttp = tps(1);
lasttp = tps(end);

H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
cellpath_name = ['/field' num2str(field) '/cellpath'];

if ~exist(fullfile(SourceF,H5filename),'file') 
    display([H5filename '.mat does not exist.']);
    return
else
    fileattrib(fullfile(SourceF,H5filename),'+w');
end

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
    H5F.close(fid);
    cellpathinfo = h5info(fullfile(SourceF,H5filename), cellpath_name);
    cellpath_mat = h5read(fullfile(SourceF,H5filename),cellpath_name,[1 1 1], [cellpathinfo.Dataspace.Size(1) cellpathinfo.Dataspace.Size(2) cellpathinfo.Dataspace.Size(3)]);
    
    for tp=firsttp:lasttp
        if ~isempty(find(cellpath_mat(:,:,tp) > 0 ,1))
            cellpath{tp} = cellpath_mat(:,:,tp);
        end
    end
else
    cellpath = [];
end
clear cellpath_mat;
datasetname = ['/field' num2str(field) '/segmentsCH' num2str(templateCH)];
selectedcells_name = ['/field' num2str(field) '/selectedcells'];
fileattrib(fullfile(SourceF,H5filename),'+w');

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,datasetname,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' datasetname]);
else
    H5L.delete(fid,datasetname,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' datasetname]);
    H5F.close(fid);
end
h5create(fullfile(SourceF,H5filename), datasetname, [size(cellpath{lasttp},1), lasttp, 3, 2*cellsize+1, 2*cellsize+1], 'Datatype', 'uint8', 'ChunkSize', [1, 1, 3, 2*cellsize+1, 2*cellsize+1], 'Deflate', 9);

selected_cells = [];
templateinfo = loadimageinfo(filetype,fileformat,[row col field plane templateCH],firsttp,channelnames,SourceF);
imwidth = templateinfo.Width;
imheight = templateinfo.Height;

%CellInd = find(cellpath{lasttp}(:,1)~=-1 & cellpath{lasttp}(:,2)~=-1)';

for cellNo = 1:size(cellpath{lasttp},1)
    for t = firsttp:lasttp
        c_cellpath(t,:) = cellpath{t}(cellNo,:);
    end
    clc;display(['r' num2str(row) 'c' num2str(col) 'f' num2str(field) ',Processing cell:' num2str(cellNo)]);
    newsegments = zeros(firsttp,lasttp,3,2*cellsize+1,2*cellsize+1, 'uint8');
    posFrameInd = 0;
    for tp = firsttp:lasttp
        
        if c_cellpath(tp,1)>0 && c_cellpath(tp,2)>0 && c_cellpath(tp,1)<= imwidth && c_cellpath(tp,2) <= imheight
            xL=max(c_cellpath(tp,1)-cellsize,1);
            xR=min(c_cellpath(tp,1)+cellsize,imwidth);
            yL=max(c_cellpath(tp,2)-cellsize,1);
            yR=min(c_cellpath(tp,2)+cellsize,imheight);
            imlocation = [row col field plane templateCH -1];
            
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
            
            template = zeros(2*cellsize+1,2*cellsize+1);
            template(borderY,borderX) = loadimage(filetype,fileformat,imlocation,tp,channelnames,{[yL yR], [xL xR]},[],SourceF);
            
            [~,~,BW] = templateToCentroid(template);
            %[~,~,BWcontour] = templateToCentroid3(template);
            if ~isempty(find(BW==1,1))
                posFrameInd = posFrameInd+1;
                newsegments(1,tp,1,:,:) = BW; %nuclei
                newsegments(1,tp,2,:,:) = bwmorph(BW,'dilate',celldilatesize); %cell
                newsegments(1,tp,3,:,:) = bwmorph(squeeze(newsegments(1,tp,2,:,:)),'erode',1) & ~squeeze(newsegments(1,tp,1,:,:)); % cytosol
            end
        end
    end
    if posFrameInd > 0.2*length(cellpath);
        selected_cells = [selected_cells;cellNo];
    end
    h5write(fullfile(SourceF,H5filename), datasetname, newsegments, [cellNo firsttp 1 1 1], [1 lasttp-firsttp+1 3 2*cellsize+1 2*cellsize+1]);
end

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,selectedcells_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' selectedcells_name]);
else
    H5L.delete(fid,selectedcells_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' selectedcells_name]);
    H5F.close(fid);
end

h5create(fullfile(SourceF,H5filename), selectedcells_name, length(selected_cells), 'Datatype', 'uint32');
h5write(fullfile(SourceF,H5filename), selectedcells_name, uint32(selected_cells));

function [x,y,BW] = templateToCentroid(M)
BWc = zeros(size(M));
for i=1.2:0.6:4
    edgedIm = edge(M,'canny',0,i);
    BW = imfill(edgedIm,'holes');
    BW = bwmorph(BW,'open',1);
    BW = bwselect(BW,round(size(M,2)/2),round(size(M,1)/2));
    BWc = BWc | BW;
end
%figure(1),imshow(BWc)
S  = regionprops(BWc, 'centroid');
BW = zeros(size(M));
if isempty(find(BWc==0, 1)) || isempty(find(BWc==1, 1))
    x = round(size(M,2)/2);
    y = round(size(M,1)/2);
else
    x = round(S.Centroid(1));
    y = round(S.Centroid(2));
    [Ys,Xs] = find(BWc==1);
    xshift = x-round(size(M,2)/2);
    yshift = y-round(size(M,1)/2);
    if isempty(find( (Ys-yshift)<=0 | (Ys-yshift)>size(M,1) | (Xs-xshift)<=0 | (Xs-xshift)>size(M,2),1))
        for i=1:length(Ys)
            BW(Ys(i)-yshift,Xs(i)-xshift) = 1;
        end
    else
        BW = BWc;
    end
end


function [x,y,BW] = templateToCentroid3(M)
xg = round(size(M,2)/2);
yg = round(size(M,1)/2);
M = imadjust(M);
M(1,:) = 0;
M(end,:) = 0;
M(:,1) = 0;
M(:,end) = 0;
outBW = modchenvese(M,80,0.05,round(size(M,1)*1));
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
boxsize = round(size(I,1)*0.1);%round(outerbox/2*.9);
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


function [notp,stagePos,stageName,waveName] = readndfile(sourcefolder,filename)
% Search for number of string matches per line.
notp=-1;
stagePos = [];
stageName = [];
waveName = [];


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


function [pathLogic,new_cellpath] = check_path(cellpath,sisterList,cellNo,firsttp,lasttp,imheight,imwidth)
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
    temp_cellpath = new_cellpath;
else
    temp_cellpath = new_cellpath(1:(deathInd-1),:);
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

function outputim = loadimage(filetype,fileformat,imlocation,tp,channelnames,pixelsreg,displayGate,SourceF)
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
            outputim = imread(fullfile(SourceF,filename),'PixelRegion',pixelsreg); 
        case 2
            outputim = imread(fullfile(SourceF,fileformat),'Index',totalCH*(tp-1)+channel1,'PixelRegion',pixelsreg);
        case 3
            filename = sprintf(fileformat,channelnames{channel1},tp);
            outputim = imread(fullfile(SourceF,filename),'PixelRegion',pixelsreg);
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
            nominim = imread(fullfile(SourceF,filename),'PixelRegion',pixelsreg);
            filename = sprintf(fileformat,row,col,field,plane,channel2,tp);
            denomin = imread(fullfile(SourceF,filename),'PixelRegion',pixelsreg);
        case 2
            nominim = imread(fullfile(SourceF,fileformat),'Index',totalCH*(tp-1)+channel1,'PixelRegion',pixelsreg);
            denomin = imread(fullfile(SourceF,fileformat),'Index',totalCH*(tp-1)+channel2,'PixelRegion',pixelsreg);
        case 3
            filename = sprintf(fileformat,channelnames{channel1},tp);
            nominim = imread(fullfile(SourceF,filename),'PixelRegion',pixelsreg);
            filename = sprintf(fileformat,channelnames{channel2},tp);
            denomin = imread(fullfile(SourceF,filename),'PixelRegion',pixelsreg);
    end
    if isempty(displayGate)
        outputim = mat2gray(im2double(nominim+2^9)./im2double(denomin+2^9));
    else
        ThresL = displayGate(1);
        ThresH = displayGate(2);
        outputim = mat2gray(im2double(nominim+2^9)./im2double(denomin+2^9),[ThresL ThresH]);
    end
end

function iminfo = loadimageinfo(filetype,fileformat,imlocation,tp,channelnames,SourceF)
row = imlocation(1);
col = imlocation(2);
field = imlocation(3);
plane = imlocation(4);
channel = imlocation(5);
totalCH = length(channelnames);

switch filetype
    case 1
        filename = sprintf(fileformat,row,col,field,plane,channel,tp);
        iminfo = imfinfo(fullfile(SourceF,filename));
    case 2
        iminfo = imfinfo(fullfile(SourceF,fileformat),'Index',totalCH*(tp-1)+channel);
    case 3
        filename = sprintf(fileformat,channelnames{channel},tp);
        iminfo = imfinfo(fullfile(SourceF,filename));
end


