function do_lsf_collectDataND130722()

% Define parameters related to the process---------
analysisparam.nucCH = 2;
analysisparam.cellCH = 2;
analysisparam.CH1 = 2;
analysisparam.CH2 = 1;
analysisparam.nominCH = 1;
analysisparam.denominCH = 2;

analysisparam.useblank_LOG = 0; % 1 = use BLANK for illumination correction, 0 = not using
analysisparam.bgnomin_LOG = 0 ; % 1 = median of BG points in SIGNAL, 2=median of BLANK
analysisparam.bgdenomin_LOG = 0; % 1 = median of BG points in SIGNAL, 2=median of BLANK
analysisparam.illumlogic = 0; % 1 = use high-pass filtering, 0 = no filtering

analysisparam.cytosize = 5;
analysisparam.cellsize = 40;
analysisparam.filterParam1 = 2;
analysisparam.filterParam2 = 2;
analysisparam.bgsize = 30;
analysisparam.signalShiftN = 0.001;
analysisparam.signalShiftD = 0.001;
analysisparam.regiontype = [2 3 2 3]; %1=nuc/cyto,2=nuc,3=cyto,4=cell
analysisparam.signaltype = [1 1 2 2];%1=ch1,2=ch2,3=nominCH/denominCH
analysisparam.regions = {'Nuc/Cyto';'Nuclei';'Cytosol';'Cell'};
analysisparam.signals = {'mCherry';'FOXO3a';'EKAREV'};
analysisparam.outputsignalNo = 1; % 'outputsignal' + 'digit'
analysisparam.output_name = 'mCherry/FOXO3a';
analysisparam.nucFolder = 'nuclearMask';
analysisparam.cellFolder = 'cellMask';
analysisparam.minstamp = 5;
analysisparam.secstamp = 0;
image_height = 1344;
image_width = 1024;

load fftexecutiontimes
h_gaussian = 1/273*[1 4 7 4 1;4 16 26 26 4;7 26 41 26 7;4 16 26 16 4;1 4 7 4 1];
smooth_opt = detbestlength2(FFTrv,FFTiv,IFFTiv,[image_height image_width],size(h_gaussian),1,1);
analysisparam.h_gaussian = h_gaussian;
analysisparam.smooth_opt = smooth_opt;
%-------------------------------------------------
% Define information about input images-----------
ndfilename ='130722.nd';
ndpathname = 'Q:\sorger\data\NIC\Bernhard\130722';
%------------------------------------------------
prefix = ndfilename(1:(end-3));
[notp,stagePos,stageName,channelnames] = readndfile(ndpathname,ndfilename);

tps = [1 notp];
sites = [64];%1:length(stagePos);

analysisparam.ndpathname = ndpathname;
analysisparam.ndfilename = ndfilename;
analysisparam.channelnames = channelnames;
analysisparam.totalCHs = 2;
analysisparam.filetype = 3;

% jobmgr = findResource('scheduler', 'type', 'lsf');
% jobmgr.ClusterMatlabRoot = '/opt/matlab';
% jobmgr.SubmitArguments = '-q short -W 12:00 -n 1 -R "rusage[matlab_dc_lic=1]" ';
% job = jobmgr.createJob();

for site = sites

    signalformat = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
    blankformat = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
    
    analysisparam.tps = tps;
    analysisparam.site = site;
    analysisparam.signalformat = signalformat;
    analysisparam.blankformat = blankformat;
    
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
    
    analysisparam.row = row;
    analysisparam.col = col;
    analysisparam.field = field;
    analysisparam.plane = plane;

    CollectData_commandline(analysisparam);
    save(fullfile(ndpathname,['site' num2str(site) 'analysis.mat']),'analysisparam');
end


function [notp,stagePos,stageName,waveName] = readndfile(ndpathname,filename)
% Search for number of string matches per line.
notp=-1;
stagePos = [];
stageName = [];
waveName = [];


if exist(fullfile(ndpathname,filename),'file')
    fid = fopen(fullfile(ndpathname,filename));
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

