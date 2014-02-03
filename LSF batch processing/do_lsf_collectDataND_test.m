function do_lsf_collectDataND_test()

% Define parameters related to the process---------
nucCH = 1;
cellCH = 2;
CH1 = 1;
CH2 = 2;
outputname = 'signals';
cytosize = 4;
cellsize = 40;
filterParam1 = 2;
filterParam2 = 2;
bgsize = 30;
signalShiftN = 0.001;
signalShiftD = 0.001;
regiontype = [1 2 3 4]; %1=nuc/cyto,2=nuc,3=cyto,4=cell
signaltype = [2 2 2 2];%4=math,1=ch1,2=ch2
regions = {'Nuc/Cyto';'Nuclei';'Cytosol';'Cell'};
signals = {'CH1';'CH2';'CH3';'Math'};
output_name = 'FOXO3a_ratio';
nucFolder = 'nuclearMask';
cellFolder = 'cellMask';
minstamp = 5;
secstamp = 0;
useblank_LOG = 0; % 1 = use BLANK for illumination correction, 0 = not using
bgnomin_LOG = 0 ; % 1 = median of BG points in SIGNAL, 2=median of BLANK
bgdenomin_LOG = 0; % 1 = median of BG points in SIGNAL, 2=median of BLANK
illumlogic = 0; % 1 = use high-pass filtering, 0 = no filtering
image_width = 1344;
image_height = 1024;

save trackingparameters;
%-------------------------------------------------
% Define information about input images-----------
ndfilename ='01262014-r1.nd';
sourcefolder = '/hms/scratch1/ss240/01-26-2014';
%------------------------------------------------
prefix = ndfilename(1:(end-3));
[notp,stagePos,stageName,channelnames] = readndfile(sourcefolder,ndfilename);
tps = [1 notp];
sites = 5;

if matlabpool('size') == 0
  matlabpool open;
end


parfor s = 1:length(sites);
    site = sites(s);
    signalformat = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
    blankformat = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
    %tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)f(?<field>\d+)','tokens');
    L = regexp(stageName{site}, 'r(?<row>\d+)','names');
    if ~isempty(L)
        row = L.row;
    else
        row = site;
    end
    
    L = regexp(stageName{site}, 'c(?<col>\d+)','names');
    if ~isempty(L)
        col = L.col;
    else
        col = 1;
    end
    L = regexp(stageName{site}, 'f(?<field>\d+)','names');
    
    if ~isempty(L)
        field = L.field;
    else
        field = 1;
    end
    
    CollectData_commandline(3,sourcefolder,row,col,field,tps,signalformat,blankformat,channelnames,2);

end


function [notp stagePos stageName waveName] = readndfile(filename)
% Search for number of string matches per line.
notp=-1;
stagePos = [];
stageName = [];
waveName = [];


if exist(filename,'file')
    fid = fopen(filename);
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


