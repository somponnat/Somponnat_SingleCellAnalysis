function do_lsf_celltrackingND_local()

% Define information about input images and necessary parameters-----------
analysisparameter.channel = 1;
analysisparameter.increment = -1;
analysisparameter.cellsize = 20;
analysisparameter.outersize = 40;
analysisparameter.similarityThres = 0.9;
maxWholeImShift = 150;
analysisparameter.maxWholeImShift = maxWholeImShift;
analysisparameter.maxNucMaskShift = 7;
analysisparameter.nucleiOptimizeLog = 1;
analysisparameter.avgNucDiameter = 17;
analysisparameter.thresParam = 8; % The higher the more stringent the threshold
analysisparameter.minAreaRatio=1.25;
analysisparameter.minCytosolWidth=5;

%---------------------------------------
ndfilename ='02272014-r1.nd';
SourceF = 'Q:\sorger\data\NIC\Pat\02-27-2014';
%------------------------------------------------
prefix = ndfilename(1:(end-3));
[notp,stagePos,stageName,channelnames] = readndfile(SourceF,ndfilename);
analysisparameter.channelnames=channelnames;
analysisparameter.SourceF=SourceF;
analysisparameter.ndfilename=ndfilename;
analysisparameter.stageName=stageName;
analysisparameter.filetype=3;

sites = 2;
tps = (notp-2):notp;
analysisparameter.tps = tps;

for i = 1:length(sites)
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
    analysisparameter.site = site;
    analysisparameter.row=row;
    analysisparameter.col=col;
    analysisparameter.field=field;
    analysisparameter.plane=plane;
    analysisparameter.fileformat=fileformat;
    
    CellTracking_commandline(analysisparameter);
    save(fullfile(SourceF,['site' num2str(site) 'tracking.mat']),'analysisparameter');
end




function [notp,stagePos,stageName,waveName] = readndfile(SourceF,filename)
% Search for number of string matches per line.
notp=-1;
stagePos = [];
stageName = [];
waveName = [];


if exist(fullfile(SourceF,filename),'file')
    fid = fopen(fullfile(SourceF,filename));
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

