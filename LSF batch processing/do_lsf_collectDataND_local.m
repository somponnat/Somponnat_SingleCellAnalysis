function do_lsf_collectDataND_local()

% Define parameters related to the process---------
analysisparam.nucCH = 1;
analysisparam.cellCH = 1;
analysisparam.CH1 = 1;
analysisparam.CH2 = 2;
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
ndfilename ='03302014-r1.nd';
ndpathname = 'Q:\sorger\data\NIC\Pat\03-30-2014';
%------------------------------------------------
prefix = ndfilename(1:(end-3));
[notp,stagePos,stageName,channelnames] = readndfile(ndpathname,ndfilename);

tps = [1 10];
sites = 1;

analysisparam.ndpathname = ndpathname;
analysisparam.ndfilename = ndfilename;
analysisparam.channelnames = channelnames;
analysisparam.totalCHs = 4;
analysisparam.filetype = 3;

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



