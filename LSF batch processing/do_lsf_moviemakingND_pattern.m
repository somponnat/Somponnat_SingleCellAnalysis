function do_lsf_moviemakingND_pattern()


% Define parameters related to the process---------
clear all;
signalshift = 0.01;
bgsubstractlogic = 0; % 
illumcorlogic = 1; % Algorithmic illumination correction by high-pass filter
framshift_logic = 0;
ImageIndex = 1; % 1=nomin/denomin, 2=templateCH, 3=nomin,4=denomin
intensityrange = [0.0032654 0.004715]; % For other images
displaygate = [0.7 1.358]; % For FRET Only
filterParam = [2 2]; 
cellsize = 15;
timestamplogic = 2; % 1 = frame no, 2 = actual time
celllocationlogic = 0; % 1 = show location of tracked cells, 0 = only image
save videoparameters;
clear all;
%-------------------------------------------------
% Define information about input images-----------
ndfilename = '05212013-r1.nd';
templateCH = 4;
nominCH = 2;
denominCH = 3;
sourcefolder = '/files/ImStor/sorger/data/NIC/Pat/05-21-2013';
%------------------------------------------------

prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(sourcefolder,ndfilename);

tps = [1 notp];
sites = 1:length(stagePos);

jobmgr = findResource('scheduler', 'type', 'lsf');
jobmgr.ClusterMatlabRoot = '/opt/matlab';
jobmgr.SubmitArguments = '-q short -W 12:00 -R "rusage[matlab_dc_lic=1]"';
job = jobmgr.createJob();

for site = sites
    
    fileformat = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
    %tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)f(?<field>\d+)','tokens');
    tokens   = regexp(stageName{site}, 'Pos(?<row>\d)-(?<col>\d+)','tokens');
    
    if ~isempty(tokens)
        row = str2num(tokens{1}{1});
        col = str2num(tokens{1}{2});
    else
        row = site;
        col = 1;
    end
    
    %field = tokens{1}{3};
    field = 1;
    plane = 1;
    job.createTask(@GenMov_commandline, 0, ...
        {3,sourcefolder, row, col,field,plane,templateCH,nominCH,denominCH, tps,fileformat,channelnames});
end

job.submit();


function [notp stagePos stageName waveName] = readndfile(sourcefolder,filename)
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
            waveName{wind} = ['w' num2str(wind) wavename1{1} '-' wavename2{1}];
            wind=wind+1;
        end
        
        tline = fgetl(fid);
    end
    fclose(fid);
end


