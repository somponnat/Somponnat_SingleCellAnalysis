function do_lsf_montagemoviemakingND_MPEG()


% Define parameters related to the process---------

signalshift = 0.001;
bgsubstractlogic = 0; % 
illumcorlogic = 1; % Algorithmic illumination correction by high-pass filter
framshift_logic = 0;
ImageIndex = 1; % 1=nomin/denomin, 2=templateCH, 3=nomin,4=denomin
intensityrange = [0.0031 0.0041]; % For other images
displaygate = [0.8 1.3]; % For FRET Only
filterParam = [2 2]; 
montageboxsize = 200;
precropscaling =0.5;
postcropscaling = 0.8;
textscaling = 1.5;
bgboxsize = 15;
videoQuality = 95;
videoFrameRate = 20;
frameoffset = [0 0];
timestamplogic = 2; % 1 = frame no, 2 = actual time
celllocationlogic = 0; % 1 = show location of tracked cells, 0 = only image
montage_name = 'EKAREV-dualsensor';
save montageparameters;

%-------------------------------------------------
% Define information about input images-----------
ndfilename = '01192014-r1.nd';

templateCH = 1;
nominCH = 2;
denominCH = 3;
sourcefolder = 'Q:\sorger\data\NIC\Pat\01-19-2014';
%------------------------------------------------

prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(sourcefolder,ndfilename);

tps = [1 150];
sites  = [16 15 14 33 34 35]; 
cellNo = [];
montageDimension = [2 3];
mInd = 1;
row = [];
col = [];
field = [];
plane = [];

for site = sites
    
    fileformat{mInd} = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
    L = regexp(stageName{site}, 'r(?<row>\d+)','names');
    if ~isempty(L)
        row(mInd) = str2num(L.row);
    else
        row(mInd) = site;
    end
    L = regexp(stageName{site}, 'c(?<col>\d+)','names');
    if ~isempty(L)
        col(mInd) = str2num(L.col);
    else
        col(mInd) = 1;
    end
    L = regexp(stageName{site}, 'f(?<field>\d+)','names');
    if ~isempty(L)
        field(mInd) = str2num(L.field);
    else
        field(mInd) = 1;
    end
    plane(mInd) = 1;
    mInd = mInd+1;
end

GenMontage_commandline(3,sourcefolder, row, col,field,plane,cellNo,templateCH,nominCH,denominCH, tps,fileformat,channelnames,montageDimension);

function [notp stagePos stageName waveName] = readndfile(sourcefolder,filename)
% Search for number of string matches per line.
notp=-1;
stagePos = [];
stageName = [];
waveName = [];

fullfile(sourcefolder,filename)

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
            wavename1  = regexp(tline, '(?<="WaveName\d+", ")\w+(?=_)', 'match')
            wavename2  = regexp(tline, '(?<="WaveName\d+", "\w+_).+(?=")', 'match')
            waveName{wind} = ['w' num2str(wind) wavename1{1} '-' wavename2{1}];
            wind=wind+1;
        end
        
        tline = fgetl(fid);
    end
    fclose(fid);
end

