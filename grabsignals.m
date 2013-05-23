function grabsignals()
clc;
close all;
ndfilename = '130220.nd';
ndpathname = 'Z:\sorger\data\NIC\Bernhard\130220_IGFRi_EGFRi\130220_IGFRi_EGFRi_Real\Convert_Files\2\';
%sourcefolder = '~/files/ImStor/sorger/data/NIC/Pat/02-03-2013/'

prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(ndpathname,ndfilename);
outputsignalNo = 1;
sequenceNo = 4;

site = 1;
tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
if ~isempty(tokens)
    row = str2num(tokens{1}{1});
    col = str2num(tokens{1}{2});
else
    row = site;
    col = 1;
end
field = 1;

H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
signal_name = ['/field' num2str(field)  '/outputsignal' num2str(outputsignalNo)];
timestamp_name = ['/field' num2str(field) '/timestamp' num2str(outputsignalNo)];
selectedcells_name = ['/field' num2str(field) '/selectedcells'];

warning off;

if exist(fullfile(ndpathname,H5filename),'file')
    fileattrib(fullfile(ndpathname,H5filename),'+w');
    fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if H5L.exists(fid,signal_name,'H5P_DEFAULT') && ...
       H5L.exists(fid,timestamp_name,'H5P_DEFAULT') 
        H5F.close(fid);
        signalinfo = h5info(fullfile(ndpathname,H5filename), signal_name);
        startind = double([1 1 sequenceNo]);
        countind = [signalinfo.Dataspace.Size(1) signalinfo.Dataspace.Size(2) 1];
        signal = permute(h5read(fullfile(ndpathname,H5filename),signal_name,startind, countind),[2 1 3]);
        timestamp = h5read(fullfile(ndpathname,H5filename),timestamp_name);

        selected_cells = 160:10:250
        for scell = 1:length(selected_cells);
            
            PosTime = find(signal(1:end,selected_cells(scell)));
            temp_y = interp1(timestamp(PosTime),signal(PosTime,selected_cells(scell)),0:5:1000);
            c_signal(:,:,scell) =  [(0:5:1000)' (1./temp_y)'];
            figure(1);
            plot(c_signal(:,1,scell),c_signal(:,2,scell),'color',[0.7 0.7 0.7]); hold on
        end
        plot(mean(c_signal(:,1,:),3),mean(c_signal(:,2,:),3),'color','k','LineWidth',2); hold on
        save(['r' num2str(row) '_c' num2str(col) '_signals'],'c_signal','selected_cells');
        assignin('base','c_signal',c_signal);
        assignin('base','selected_cells',selected_cells);
    end
end



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

