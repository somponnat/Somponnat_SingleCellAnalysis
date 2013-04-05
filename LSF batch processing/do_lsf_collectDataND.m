function do_lsf_collectDataND()

% Define parameters related to the process---------
clear all;
CH1 = 4;
CH2 = 2;
CH3 = 3;
nominCH = 2;
denominCH = 3;
outputname = 'signals';
cytosize = 4;
cellsize = 40;
filterParam1 = 2;
filterParam2 = 2;
bgsize = 30;
signalShiftN = 0.001;
signalShiftD = 0.001;
regiontype = [2 3 2 3]; %1=nuc/cyto,2=nuc,3=cyto,4=cell
signaltype = [4 4 1 1];%4=math,1=ch1,2=ch2,3=ch3
minstamp = 5;
secstamp = 0;
useblank_LOG = 1; % 1 = use BLANK for illumination correction, 0 = not using
bgnomin_LOG = 2 ; % 1 = median of BG points in SIGNAL, 2=median of BLANK
bgdenomin_LOG = 2; % 1 = median of BG points in SIGNAL, 2=median of BLANK
illumlogic = 0; % 1 = use high-pass filtering, 0 = no filtering
image_width = 1344;
image_height = 1024;

save trackingparameters;
clear all;
%-------------------------------------------------
% Define information about input images-----------
ndfilename = '02032013-r1.nd';
templateCH = 1;
BLANKsite = 25;
sourcefolder = '/files/ImStor/sorger/data/NIC/Pat/02-03-2013';
%------------------------------------------------
currentF = pwd;
cd(sourcefolder);
prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(ndfilename);
cd(currentF);
tps = [1 notp];
sites = 1:length(stagePos);

jobmgr = findResource('scheduler', 'type', 'lsf');
jobmgr.ClusterMatlabRoot = '/opt/matlab';
jobmgr.SubmitArguments = '-q short -W 12:00 -R "rusage[matlab_dc_lic=1]"';
job = jobmgr.createJob();

for site = sites
    
    signalformat = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
    blankformat = [prefix '_%s_s' num2str(BLANKsite) '_t%g.TIF'];
    
    %tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)f(?<field>\d+)','tokens');
    tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
    
    if ~isempty(tokens)
        row = str2num(tokens{1}{1});
        col = str2num(tokens{1}{2});
    else
        row = site;
        col = 1;
    end
    
    %field = tokens{1}{3};
    field = 1;
    
    job.createTask(@CollectData_commandline, 0, ...
        {3,sourcefolder,row,col,field,tps,signalformat,blankformat,channelnames,totalCHs,templateCH});

end

job.submit();

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


