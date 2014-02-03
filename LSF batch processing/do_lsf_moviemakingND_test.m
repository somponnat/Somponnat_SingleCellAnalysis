function do_lsf_moviemakingND_test()


% Define parameters related to the process---------
clear all;
signalshift = 0.001;
bgsubstractlogic = 0; % 
illumcorlogic = 1; % Algorithmic illumination correction by high-pass filter
framshift_logic = 0;
ImageIndex = 2; % 1=nomin/denomin, 2=templateCH, 3=nomin,4=denomin
intensityrange = [0.0031281 0.0042]; % For other images
displaygate = [0.27 2.316]; % For FRET Only
filterParam = [2 2]; 
cellsize = 15;
timestamplogic = 2; % 1 = frame no, 2 = actual time
celllocationlogic = 0; % 1 = show location of tracked cells, 0 = only image
save videoparameters;
clear all;
%-------------------------------------------------
% Define information about input images-----------
ndfilename = '01182014-r2.nd';
templateCH = 1;
nominCH = 2;
denominCH = 3;
sourcefolder = 'Q:\sorger\data\NIC\Pat\01-18-2014';
%------------------------------------------------

prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(sourcefolder,ndfilename);

tps = [1 notp];
sites = 1:length(stagePos);

if matlabpool('size') == 0
  matlabpool open;
end

parfor site = sites
    
    fileformat = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
    tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)f(?<field>\d+)','tokens');
    row = tokens{1}{1};
    col = tokens{1}{2};
    field = tokens{1}{3};
    plane = 1;
    GenMov_commandline(3,sourcefolder, row, col,field,plane,templateCH,nominCH,denominCH, tps,fileformat,channelnames);
end

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

