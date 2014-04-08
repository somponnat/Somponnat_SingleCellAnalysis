function do_lsf_collectDataND()

% Define parameters related to the process---------
nucCH = 1;
cellCH = 2;
CH1 = 1;
CH2 = 2;
nominCH = 1;
denominCH = 2;

useblank_LOG = 0; % 1 = use BLANK for illumination correction, 0 = not using
bgnomin_LOG = 0 ; % 1 = median of BG points in SIGNAL, 2=median of BLANK
bgdenomin_LOG = 0; % 1 = median of BG points in SIGNAL, 2=median of BLANK
illumlogic = 0; % 1 = use high-pass filtering, 0 = no filtering

cytosize = 5;
cellsize = 40;
filterParam1 = 2;
filterParam2 = 2;
bgsize = 30;
signalShiftN = 0.001;
signalShiftD = 0.001;
regiontype = [2 3 2 3]; %1=nuc/cyto,2=nuc,3=cyto,4=cell
signaltype = [1 1 2 2];%1=ch1,2=ch2,3=nominCH/denominCH
regions = {'Nuc/Cyto';'Nuclei';'Cytosol';'Cell'};
signals = {'mCherry';'FOXO3a';'EKAREV'};
outputsignalNo = 1; % 'outputsignal' + 'digit'
output_name = 'mCherry/FOXO3a';
nucFolder = 'nuclearMask';
cellFolder = 'cellMask';

minstamp = 5;
secstamp = 0;
image_width = 1344;
image_height = 1024;

save trackingparameters;
%-------------------------------------------------
% Define information about input images-----------
ndfilename ='02122014wtakt-r2.nd';
sourcefolder = '/hms/scratch1/ss240/02-12-2014-wtAkt';
%------------------------------------------------
prefix = ndfilename(1:(end-3));
[notp,stagePos,stageName,channelnames] = readndfile(sourcefolder,ndfilename);
tps = [1 notp];
sites = 1:length(stagePos);

jobmgr = findResource('scheduler', 'type', 'lsf');
jobmgr.ClusterMatlabRoot = '/opt/matlab';
jobmgr.SubmitArguments = '-q short -W 12:00 -n 1 -R "rusage[matlab_dc_lic=1]" ';
job = jobmgr.createJob();

for site = sites
    
    signalformat = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
    blankformat = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
    
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
    

    job.createTask(@CollectData_commandline, 0, ...
        {3,sourcefolder,row,col,field,tps,signalformat,blankformat,channelnames,2});
end

job.submit();

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



