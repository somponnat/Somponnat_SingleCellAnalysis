function do_lsf_celltrackingND3()

% Define information about input images and necessary parameters-----------
templateCH = 1;
increment = -1;
cellsize = 18;
outersize = 40;
similarityThres = 0.9;
maxWholeImShift = 100;
maxNucMaskShift = 15;
nucleiOptimizeLog = 1;
avgNucDiameter = 17;
thresParam = 8; % The higher the more stringent the threshold
minAreaRatio=1.25;
minCytosolWidth=5;
save celltrackingparameters4;
%---------------------------------------
ndfilename ='03312014-r2.nd';
sourcefolder = '/hms/scratch1/ss240/03-31-2014';
%------------------------------------------------
prefix = ndfilename(1:(end-3));
[notp,stagePos,stageName,channelnames] = readndfile(sourcefolder,ndfilename);
tps = [1 notp];
sites = [9	10	11	12 16	15	14	13 33	34	35	36 40	39	38	37 57	58	59	60 64	63	62	61];

jobmgr = findResource('scheduler', 'type', 'lsf');
jobmgr.ClusterMatlabRoot = '/opt/matlab';
jobmgr.SubmitArguments = '-q short -W 12:00 -M 16777216 -n 1 -R "rusage[matlab_dc_lic=1]"';
job = jobmgr.createJob();

for site = sites
    
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
    
    job.createTask(@CellTracking_commandline3, 0, ...
        {3,sourcefolder,row,col,field,plane,templateCH,tps,increment,fileformat,channelnames,cellsize,outersize,similarityThres});

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
