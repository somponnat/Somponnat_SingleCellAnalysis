function do_lsf_imagesegment_local()

% Define information about input images and necessary parameters-----------
segmentparameter.channel = 1;
segmentparameter.increment = 1;
segmentparameter.thresParam = 2; % The higher the more stringent the threshold
segmentparameter.minNucDiameter = 11;
segmentparameter.maxNucDiameter = 130;
segmentparameter.thresParam = 3;
segmentparameter.minAreaRatio = 1.05;
segmentparameter.minCytosolWidth = 5;
%---------------------------------------
ndfilename = '03292014-r1.nd';
SourceF = 'Q:\sorger\data\NIC\Pat\03-29-2014';
%------------------------------------------------
prefix = ndfilename(1:(end-3));
[notp,stagePos,stageName,channelnames] = readndfile(SourceF,ndfilename);
segmentparameter.channelnames=channelnames;
segmentparameter.SourceF=SourceF;
segmentparameter.ndfilename=ndfilename;
segmentparameter.stageName=stageName;
segmentparameter.filetype=3;

sites  = 1:length(stagePos);
tps = [1 notp];

for i = 1:length(sites)
    site = sites(i);
    segmentparameter.tps = tps;
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
    segmentparameter.site = site;
    segmentparameter.row=row;
    segmentparameter.col=col;
    segmentparameter.field=field;
    segmentparameter.plane=plane;
    segmentparameter.fileformat=fileformat;
    
    ImageSegmenting_commandline(segmentparameter);
    save(fullfile(SourceF,['site' num2str(site) 'segmentparam.mat']),'segmentparameter');
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
