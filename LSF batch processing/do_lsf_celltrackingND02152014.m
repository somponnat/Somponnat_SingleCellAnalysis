function do_lsf_celltrackingND02152014()

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
ndfilename ='02152014-r3.nd';
SourceF = '/hms/scratch1/ss240/02-15-2014';
%------------------------------------------------
prefix = ndfilename(1:(end-3));
[notp,stagePos,stageName,channelnames] = readndfile(SourceF,ndfilename);
analysisparameter.channelnames=channelnames;
analysisparameter.SourceF=SourceF;
analysisparameter.ndfilename=ndfilename;
analysisparameter.stageName=stageName;
analysisparameter.filetype=3;

sites = [14];

tps = [1 notp];
analysisparameter.tps = tps;

jobmgr = findResource('scheduler', 'type', 'lsf');
jobmgr.ClusterMatlabRoot = '/opt/matlab';
jobmgr.SubmitArguments = '-q sysbio_1d -n 1 -R "rusage[matlab_dc_lic=1]"';
job = jobmgr.createJob();

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
    
    job.createTask(@CellTracking_commandline, 0, {analysisparameter});
    save(fullfile(SourceF,['site' num2str(site) 'tracking.mat']),'analysisparameter');
end

job.submit();


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
