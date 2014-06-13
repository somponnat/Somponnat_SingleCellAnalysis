function ImageSegmenting_commandline(segmentparameter)

filetype = segmentparameter.filetype;
SourceF = segmentparameter.SourceF;
row = segmentparameter.row;
col = segmentparameter.col;
field = segmentparameter.field;
plane = segmentparameter.plane;
channel = segmentparameter.channel;
tps = segmentparameter.tps;
increment = segmentparameter.increment;
fileformat = segmentparameter.fileformat;
channelnames = segmentparameter.channelnames;
minNucDiameter = segmentparameter.minNucDiameter;
maxNucDiameter = segmentparameter.maxNucDiameter;
thresParam = segmentparameter.thresParam;
minAreaRatio = segmentparameter.minAreaRatio;
minCytosolWidth = segmentparameter.minCytosolWidth;

currentPath = pwd;
eval('cd ..');
addpath(genpath([pwd filesep 'ThirdParty']),'-end');
cd(currentPath);
switch increment
    case 1
        firstFrame = tps(1);
        endFrame   = tps(end);
    case -1
        firstFrame = tps(end);
        endFrame   = tps(1);
end

H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
cellpath_name = ['/field' num2str(field) '/cellpath'];


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
padWidth =15;
minFormfactor=0.7;
for tp = firstFrame:increment:endFrame
    
    currentframe = loadimage(filetype,fileformat,[row col field plane channel],tp,channelnames,SourceF);

    [ySize, xSize] = size(currentframe);
    currentframe=im2uint16(currentframe);
    [iMax, iWidth] = contrastValues(currentframe);
    tmp = iMax + (rand(ySize+2*padWidth, xSize+2*padWidth)-0.5) * 2*iWidth;
    tmp(padWidth+1:ySize+padWidth, padWidth+1:xSize+padWidth) = currentframe;
    currentframe = tmp/max(tmp(:));
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
    currentframe = currentframe(padWidth+1:ySize+padWidth, padWidth+1:xSize+padWidth);
    
    % define cell boundary
    bwCell=cellSegmentation(currentframe,bwNuc,thresParam,minAreaRatio,minCytosolWidth);
    
    % determine centroid
    CC = bwconncomp(bwNuc, 8);
    S = regionprops(CC,'Centroid');
    numCell = length(S);
    
    if numCell>0
        cellpath_mat = -1*(ones(numCell,2,tps(end)));
        nucCent=zeros(size(currentframe));
        for c=1:numCell
            cellpath_mat(c,1,tp) = round(S(c).Centroid(1));
            cellpath_mat(c,2,tp) = round(S(c).Centroid(2));
            nucCent(cellpath_mat(c,2,tp),cellpath_mat(c,1,tp))=1;
        end
        
        p2=imdilate(nucCent,strel('disk',2));
        p3=bwperim(bwNuc);
        p3=p3|bwperim(bwNuc&~p3);
        p4=bwperim(bwNuc);
        p4=p4|bwperim(bwCell&~p4);
        im=imadjust(im2double(currentframe));
        imOut=im2uint8(cat(3,max(im-p2-p3-p4,p2),max(im-p2-p3-p4,p3),max(im-p2-p3-p4,p4)));
        writeimage(imOut,filetype,fileformat,[row col field plane channel],tp,channelnames,overlayF);
        writeimage(bwNuc,filetype,fileformat,[row col field plane channel],tp,channelnames,nucMaskF);
        writeimage(bwCell,filetype,fileformat,[row col field plane channel],tp,channelnames,cellMaskF);
        
        display(['Successfully seeded frame ' num2str(tp) ' and saved in ' H5filename]);
    end
    
end

function BW = readBWimage(filetype,fileformat,imlocation,tp,channelnames,destF)
row = imlocation(1);
col = imlocation(2);
field = imlocation(3);
plane = imlocation(4);
CH = imlocation(5);

switch filetype
    case 1
        
        filename = sprintf(fileformat,row,col,field,plane,CH,tp);
        BW = imread(fullfile(destF,filename));
    case 2
        
    case 3
        filename = sprintf(fileformat,channelnames{CH},tp);
        BW = imread(fullfile(destF,filename));
        
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

function out = checkimage(filetype,fileformat,imlocation,tp,channelnames,destF)
row = imlocation(1);
col = imlocation(2);
field = imlocation(3);
plane = imlocation(4);
CH = imlocation(5);

switch filetype
    case 1
        
        filename = sprintf(fileformat,row,col,field,plane,CH,tp);
        out = exist(fullfile(destF,filename),'file');
    case 3
        filename = sprintf(fileformat,channelnames{CH},tp);
        out = exist(fullfile(destF,filename),'file');
        
end


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
bwCell=dilateNuc ;%| bwCell;
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

function [x,y,BW] = templateToCentroid(M,xg,yg)
BWc = zeros(size(M));
%oldArea = [];
for i=1.2:0.4:6
    
    edgedIm = edge(M,'canny',0,i);
    BW = imfill(edgedIm,'holes');
    BW = bwmorph(BW,'open',1);
    BW = bwselect(BW,xg,yg);
    newArea = bwarea(BW);
    BWc = BWc | BW;
    
end

BW = BWc;
S  = regionprops(BW, 'centroid');

if isempty(find(BW==0, 1)) || isempty(find(BW==1, 1))
    x = [];
    y = [];
    
else
    x = round(S.Centroid(1));
    y = round(S.Centroid(2));
end

function [x,y,BW] = templateToCentroid3(M,startBW,xg,yg)
ffactorCutoff = 0.7;
M(1,:) = 0;
M(end,:) = 0;
M(:,1) = 0;
M(:,end) = 0;
outBW = modchenvese(M,startBW,100,0.08,28);
se = strel('disk',1);
outBW = imclose(outBW,se);

if outBW(yg,xg)==0
    BW = false(size(M));
    x = xg;
    y = yg;
    return;
end
outBW(1,:) = false;
outBW(end,:) = false;
outBW(:,1) = false;
outBW(:,end) = false;
outBW = imfill(outBW,'holes');
selectedBW = bwselect(outBW,xg,yg);
selectedNEGBW = bwselect(~outBW,xg,yg);
S1  = regionprops(selectedBW, 'centroid','Area','Perimeter');
S2  = regionprops(selectedNEGBW, 'centroid','Area','Perimeter');
S3  = regionprops(outBW, 'centroid','Area','Perimeter');
S4  = regionprops(~outBW, 'centroid','Area','Perimeter');

if ~isfield(S1,'Centroid') || S1.Area > 0.8*length(M)^2
    x = xg;
    y = yg;
    BW = false(size(M));
elseif length(S1)==1
    if 4*pi*S1.Area/S1.Perimeter^2 > ffactorCutoff
        x = round(S1.Centroid(1));
        y = round(S1.Centroid(2));
        BW = selectedBW;
    else
        x = xg;
        y = yg;
        BW = false(size(M));
    end
elseif length(S2)==1
    if 4*pi*S2.Area/S2.Perimeter^2 > ffactorCutoff
        x = round(S2.Centroid(1));
        y = round(S2.Centroid(2));
        BW = selectedNEGBW;
    else
        x = xg;
        y = yg;
        BW = false(size(M));
    end
elseif length(S3)==1 && 4*pi*S3.Area/S3.Perimeter^2 > ffactorCutoff
    x = round(S3.Centroid(1));
    y = round(S3.Centroid(2));
    BW = outBW;
elseif length(S4)==1 && 4*pi*S4.Area/S4.Perimeter^2 > ffactorCutoff
    x = round(S4.Centroid(1));
    y = round(S4.Centroid(2));
    BW = ~outBW;
elseif length(S3)>1
    currentD = [];
    for i=1:length(S3)
        if 4*pi*S3(i).Area/S3(i).Perimeter^2 > ffactorCutoff
            currentD(i) = (S3(i).Centroid(1)-xg)^2 + (S3(i).Centroid(2)-yg)^2;
        else
            currentD(i) = 1000;
        end
    end
    goodInd = find( currentD == min(currentD) );
    if ~isempty(goodInd)
        x=S3(goodInd).Centroid(1);
        y=S3(goodInd).Centroid(2);
        BW = selectedBW(outBW,x,y);
    else
        x = xg;
        y = yg;
        BW = false(size(M));
    end
elseif length(S4)>1
    currentD = [];
    for i=1:length(S4)
        if 4*pi*S4(i).Area/S4(i).Perimeter^2 > ffactorCutoff
            currentD(i) = (S4(i).Centroid(1)-xg)^2 + (S4(i).Centroid(2)-yg)^2;
        else
            currentD(i) = 1000;
        end
    end
    goodInd = find( currentD == min(currentD) );
    if ~isempty(goodInd)
        x=S4(goodInd).Centroid(1);
        y=S4(goodInd).Centroid(2);
        BW = selectedBW(~outBW,x,y);
    else
        x = xg;
        y = yg;
        BW = false(size(M));
    end
else
    x = xg;
    y = yg;
    BW = false(size(M));
end

%   Adapted from code by Yue Wu (yue.wu@tufts.edu)
%   http://www.mathworks.com/matlabcentral/fileexchange/23445
function seg = modchenvese(I,startBW,num_iter,mu,outerbox)
s = outerbox./min(size(I,1),size(I,2)); % resize scale
I = imresize(I,s);
mask = imresize(startBW,s);

phi0 = bwdist(mask)-bwdist(1-mask)+im2double(mask)-.5;
force = eps;

%-- Display settings
%figure();
%subplot(2,2,1); imshow(I); title('Input Image');axis equal;
%subplot(2,2,2); contour(flipud(phi0), 'r','LineWidth',1); title('initial contour');axis equal;
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
    %   showphi(I,phi0,n);
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
%figure(3);imshow(corrScore); drawnow;
% 3. finding most likely region

[maxVal,maxIdx] = max(corrScore(:));
[y, x] = ind2sub([size(corrScore,1),size(corrScore,2)],maxIdx);
%figure(3); hold on; plot(x,y,'rx');
x = x-round(templateWidth/2);
y = y-round(templateHeight/2);

