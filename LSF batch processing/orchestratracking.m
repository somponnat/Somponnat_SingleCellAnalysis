function orchestratracking()

rows = 2:3;
cols = 2:4;
fields = 1:4;
channels = 4;
plane = 1;
tps = [1 130];
lind = 1;

matlabpool open;

for row = rows
    for col = cols
        for field = fields
            for channel = channels
                loopparam(lind,:) = [row col field channel];
                lind = lind+1;
            end
        end
    end
end


parfor lind = 1:size(loopparam,1)
    row = loopparam(lind,1);
    col = loopparam(lind,2);
    field = loopparam(lind,3);
    channel = loopparam(lind,4);
    CellTracking_commandline(row,col,field,plane,channel,tps);
end
matlabpool close;

function CellTracking_commandline(row,col,field,plane,channel,tps)
% Example usage: CellTracking_commandline(pwd,2,6,1,1,3,[1 16]);
% celltrackOUT    =  MAT file that contains initial track points
% tps             = [<first frame>   <last frame>]
%------------------------------
% Adjustable parameters
outersize = 65;
cellsize = 40;
fileformat = 'r%02.0fc%02.0ff%02.0fp%02.0frc%1.0f-ch1sk%ufk1fl1.tiff';
similarityThres = 0.8;
maxI = 1e11;
invertedLog = 0;
cropFactor = 0.1;
%------------------------------
currentF = pwd;
load fftexecutiontimes
warning('off', 'all');

savefile = ['celltrackOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(channel)];

% check if necessary file exist, if not quit jobs
if ~exist([savefile '.mat'],'file')
    display([savefile '.mat does not exist.']);
    return
end
% file exist, load the file
load(savefile,'cellpath','sisterList','bg','mark');

% check if the necessary variables exist, if not quit jobs
if ~exist('cellpath','var') | ~exist('bg','var') | ~exist('sisterList','var')
    display([savefile '.mat does not contain necessary variables.']);
    return
end

% calculate optimal conditions for fft analysis
fftw('planner', 'hybrid');
opt = detbestlength2(FFTrv,FFTiv,IFFTiv,2*[outersize outersize],2*[cellsize cellsize],1,1);

for t=tps(1):tps(2)-1
    clc;display([savefile ' - Currently processing frame: ' num2str(t+1) ' of ' num2str(tps(2)-tps(1)+1)]);
    tp=t; % load previous file
    filename = sprintf(fileformat,row,col,field,plane,channel,tp);
    previousframe = imread(filename);

    tp=t+1; % load current file
    filename = sprintf(fileformat,row,col,field,plane,channel,tp);
    currentframe = imread(filename);
    
    sisterList{tp} = sisterList{tp-1};

    for cell=1:size(cellpath{tp-1},1)
        
        xL=max(cellpath{tp-1}(cell,1)-cellsize,1);
        xR=min(cellpath{tp-1}(cell,1)+cellsize-1,size(previousframe,2));
        yL=max(cellpath{tp-1}(cell,2)-cellsize,1);
        yR=min(cellpath{tp-1}(cell,2)+cellsize-1,size(previousframe,1));
                
        template = previousframe(yL:yR,xL:xR);     
        
        xL=max(cellpath{tp-1}(cell,1)-outersize,1);
        xR=min(cellpath{tp-1}(cell,1)+outersize-1,size(currentframe,2));
        yL=max(cellpath{tp-1}(cell,2)-outersize,1);
        yR=min(cellpath{tp-1}(cell,2)+outersize-1,size(currentframe,1));
        
        testframe = currentframe(yL:yR,xL:xR);   
        
        [x1 y1] = corrMatching2(testframe, template,similarityThres,opt);
        [x1 y1 BW] = templateToCentroid(testframe,x1,y1,maxI,invertedLog);

        cellpath{tp}(cell,:) = [xL+x1 yL+y1]; 
    end
    
    if ~isempty(bg)
        %updating bg points
        xshift = mean(cellpath{tp}(:,1))-mean(cellpath{tp-1}(:,1));
        yshift = mean(cellpath{tp}(:,2))-mean(cellpath{tp-1}(:,2));
        
        bg{tp}(:,1) = bg{tp-1}(:,1)+round(xshift);
        bg{tp}(:,2) = bg{tp-1}(:,2)+round(yshift);

    end

end

% finishing jobs, returning to original path
save(savefile,'cellpath','sisterList','bg','mark','-v7.3');
display(['Successfully tracked frame' num2str(tps(1)) ':' num2str(tps(2)) ' and saved in ' savefile]);
cd(currentF);

function [x y BW] = templateToCentroid(M,xg,yg,maxI,invertLog)
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



function [x y] = corrMatching2(frameImg, templateImg, threshC,opt)
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

