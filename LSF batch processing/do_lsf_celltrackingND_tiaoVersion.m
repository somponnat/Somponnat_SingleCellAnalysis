function do_lsf_celltrackingND_tiaoVersion()
% clc;
% clear all;
% Define information about input images and necessary parameters-----------
cellsize = 15;
outersize = 60;
avgNucDiameter = 12;
maxWholeImShift = 300;
maxNucMaskShift = 10;
similarityThres = 0.9;
ffactorCutoff = 0.85; %roundness
thresParam = 3; % The higher the more stringent the threshold
minAreaRatio=1.25;
minCytosolWidth=5;
save celltrackingparameters;

%-------------------------------------------------
% Define information about input images-----------
ndfilename ='02272014-r1.nd';
sourcefolder = '/hms/scratch1/ss240/02-27-2014';
%------------------------------------------------
prefix = ndfilename(1:(end-3));
[notp,stagePos,stageName,channelnames] = readndfile(sourcefolder,ndfilename);
tps = [1 notp];
%sites = 1:length(stagePos);
sites = [1	2	3	4	5	6,...
18	17	16	15	14	13,...
19	20	21	22	23	24,...
36	35	34	33	32	31,...
37	38	39	40	41	42,...
54	53	52	51	50	49,...
55	56	57	58	59	60];
jobmgr = findResource('scheduler', 'type', 'lsf');
jobmgr.ClusterMatlabRoot = '/opt/matlab';
jobmgr.SubmitArguments = '-q sysbio_7d -M 16777216 -n 1 -R "rusage[matlab_dc_lic=1]" ';
job = jobmgr.createJob();

for site = sites
    NucCH = 1;
    CellCH = 2;
    increment = 1;
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

    job.createTask(@CellTracking_commandline_tiao, 0, ...
        {3,sourcefolder,row,col,field,plane,NucCH,CellCH,tps,increment,fileformat,channelnames});
    
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



