function createClusterParams()
clc;

% Define location of HDF5 files and the original ND file
ndfilename = '02032013-r1.nd';
%sourcefolder = 'C:\computation\02-03-2013';
sourcefolder = 'Q:\sorger\data\computation\Bernhard_Steiert\EKAREV dynamics\02-03-2013-dataanalysis';

prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(sourcefolder,ndfilename); % read in information about the experiment

% the output signal dataset number
outputsignalNo = 1;

% signal type within the dataset
sequenceNo = 2;

% specific the site(s) to be processed
sites = 7;

% check with matlabpool is already initiated
% if matlabpool('size') == 0
%   matlabpool open;
% end

% loop through all sites and dispatch jobs using parfor
for s=1:length(sites)
    site = sites(s);
    tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
    if ~isempty(tokens)
        row = str2num(tokens{1}{1});
        col = str2num(tokens{1}{2});
    else
        row = site;
        col = 1;
    end
    field = 1;
    cal_clusterparam(row,col,field,sourcefolder,outputsignalNo,sequenceNo);
end


% This is the main function that preprocess signal and convert original
% time-series to parameters.
function cal_clusterparam(row,col,field,ndpathname,outputsignalNo,sequenceNo)
minimumSignalSize = 100;
MiddleToTop = 1; % 1 = assugb the middle cluster to top group, 0 = assign the middle cluster to bottom group
showPlots = 1; % change to 1 if needing to visualize the peak detection
midHgating = 0.08; % x the median of lowest peak cluster
delayGate = 0.5; % fraction of height that must decay to consider as peak tail

% specifying the pre- and post-treatment time
PreTime = 1:25;
PostTime = 40:244;

% this section specific the naming of HDF5 file as well as subdirectories
H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
signal_name = ['/field' num2str(field)  '/outputsignal' num2str(outputsignalNo)];
timestamp_name = ['/field' num2str(field) '/timestamp' num2str(outputsignalNo)];
selectedcells_name = ['/field' num2str(field) '/selectedcells'];
param_name = ['/field' num2str(field)  '/clusterparams' num2str(outputsignalNo)];
cellpath_name = ['/field' num2str(field) '/cellpath'];
sisterList_name = ['/field' num2str(field) '/sisterList'];
division_name = ['/field' num2str(field)  '/divisiontime'];
peak_name = ['/field' num2str(field)  '/peakmat' num2str(outputsignalNo)];

warning off;

% first check if the HDF5 file already exists
if exist(fullfile(ndpathname,H5filename),'file')
    % If so, enables the write ability
    fileattrib(fullfile(ndpathname,H5filename),'+w');
    % create file handler using low-level h5 function
    fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    % check with time-series, timestamp, and selected matrices exist
    if H5L.exists(fid,signal_name,'H5P_DEFAULT') && ...
            H5L.exists(fid,timestamp_name,'H5P_DEFAULT') && ...
            H5L.exists(fid,selectedcells_name,'H5P_DEFAULT')
        % need to explicitly close the file handler when using low-level
        % h5 function
        H5F.close(fid);
        
        % determine the size of the time-series matrix
        signalinfo = h5info(fullfile(ndpathname,H5filename), signal_name);
        startind = double([1 1 sequenceNo]); % note here that we only grab the specified slice of sequenceNo
        countind = [signalinfo.Dataspace.Size(1) signalinfo.Dataspace.Size(2) 1];
        
        % permute dimension of the time-series signal matrix to:
        % Timepoint,cellNo, sequenceNo
        signal = permute(h5read(fullfile(ndpathname,H5filename),signal_name,startind, countind),[2 1 3]);
        timestamp = h5read(fullfile(ndpathname,H5filename),timestamp_name);
        selected_cells = h5read(fullfile(ndpathname,H5filename),selectedcells_name);

        cellpathinfo = h5info(fullfile(ndpathname,H5filename), cellpath_name);
        sisterListinfo = h5info(fullfile(ndpathname,H5filename), sisterList_name);
        
        % read in cellpath and sisterList matrices
        cellpath_mat = h5read(fullfile(ndpathname,H5filename),cellpath_name,[1 1 1], [cellpathinfo.Dataspace.Size(1) cellpathinfo.Dataspace.Size(2) cellpathinfo.Dataspace.Size(3)]);
        sisterList_mat = h5read(fullfile(ndpathname,H5filename),sisterList_name,[1 1 1], [sisterListinfo.Dataspace.Size(1) sisterListinfo.Dataspace.Size(2) sisterListinfo.Dataspace.Size(3)]);
        
        % Initialize peaks matrix
        fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
        if ~H5L.exists(fid,peak_name,'H5P_DEFAULT')
            H5F.close(fid);
            display(['Initializing ' H5filename ':' peak_name]);
        else
            H5L.delete(fid,peak_name,'H5P_DEFAULT');
            display(['Overwriting ' H5filename ':' peak_name]);
            H5F.close(fid);
        end
        
        % cellNo, peak type, peak params, peak#
        h5create(fullfile(ndpathname,H5filename), peak_name, [size(cellpath_mat,1),2, 4, 200], 'Datatype', 'double');
        % cellNo, division time (1-3 rounds)
 
        fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
        if ~H5L.exists(fid,division_name,'H5P_DEFAULT')
            H5F.close(fid);
            display(['Initializing ' H5filename ':' division_name]);
        else
            H5L.delete(fid,division_name,'H5P_DEFAULT');
            display(['Overwriting ' H5filename ':' division_name]);
            H5F.close(fid);
        end
        h5create(fullfile(ndpathname,H5filename), division_name, [size(cellpath_mat,1),3], 'Datatype', 'double');

        param_mat = [];
        
        signal_mean_t = smooth(mean(signal,2));
        
        % loop through all cells in the selected-cell lists
        for scell = setdiff(selected_cells',[]); %%%%
            
            % determine the phenotype of the cells
            if cellpath_mat(scell,1,end) > 0
                if     sisterList_mat(scell,1,end) ~= -1 && sisterList_mat(scell,2,end) == -1 && sisterList_mat(scell,3,end) == -1
                    phenotype = 1; %dividing once
                elseif sisterList_mat(scell,1,end) ~= -1 && sisterList_mat(scell,2,end) ~= -1 && sisterList_mat(scell,3,end) == -1
                    phenotype = 2; %dividing twice
                elseif sisterList_mat(scell,1,end) ~= -1 && sisterList_mat(scell,2,end) ~= -1 && sisterList_mat(scell,3,end) ~= -1
                    phenotype = 3; %dividing three times
                else
                    phenotype = 0; %non-dividing
                end
            elseif cellpath_mat(scell,1,end) == -2
                if  sisterList_mat(scell,1,end) ~= -1 && sisterList_mat(scell,2,end) == -1 && sisterList_mat(scell,3,end) == -1
                    phenotype = -1; %dividing once
                elseif sisterList_mat(scell,1,end) ~= -1 && sisterList_mat(scell,2,end) ~= -1 && sisterList_mat(scell,3,end) == -1
                    phenotype = -2; %dividing twice
                elseif sisterList_mat(scell,1,end) ~= -1 && sisterList_mat(scell,2,end) ~= -1 && sisterList_mat(scell,3,end) ~= -1
                    phenotype = -3; %dividing three times
                else
                    phenotype = -4; %non-dividing
                end
                display('cell dead');
            end
            
            %find division time
            divisionTime = zeros(1,3,'double');
            switch phenotype
                case {1,-1}
                    divisionTime(1) = timestamp(find(sisterList_mat(scell,1,:) ~= -1,1,'first'));
                case {2,-2}
                    divisionTime(1) = timestamp(find(sisterList_mat(scell,1,:) ~= -1,1,'first'));
                    divisionTime(2) = timestamp(find(sisterList_mat(scell,2,:) ~= -1,1,'first'));
                case {3,-3}
                    divisionTime(1) = timestamp(find(sisterList_mat(scell,1,:) ~= -1,1,'first'));
                    divisionTime(2) = timestamp(find(sisterList_mat(scell,2,:) ~= -1,1,'first'));
                    divisionTime(3) = timestamp(find(sisterList_mat(scell,3,:) ~= -1,1,'first'));
            end
            
            %save information about cell division time
            h5write(fullfile(ndpathname,H5filename), division_name, divisionTime, [double(scell) 1], [1 3]);
            
            %find non-zero timepoint
            PosTime = find(signal(:,scell)~=0);
            c_signal = signal(PosTime,scell);
            c_time   = timestamp(PosTime);
            
            %only process further if the time-series have minimum signal
            %size longer than the minimum signal size
            if numel(PosTime)>minimumSignalSize
                display([H5filename 'cell:' num2str(scell)]);
                xs = c_time(1):7:timestamp(PosTime(end));
                ys_ori=interp1(c_time,c_signal,xs);
                ys=ys_ori;
                
                signal_mean=interp1(timestamp,signal_mean_t,xs);
                ys = ys_ori-signal_mean;
                outTS = getTimeSeriesTrend(ys,'trendType',1);
                ys = outTS.dTS;
                
                % determining and assigning peaks to blank 4-by-200 matrix
                p_peaks = zeros(4,200,'double');
                t_peaks = zeros(4,200,'double');
                p_temp = findTruePeaks(xs,ys,0,MiddleToTop,midHgating,delayGate,[timestamp(PreTime(1)) timestamp(PreTime(end))],divisionTime);
                t_temp = findTruePeaks(xs,ys,showPlots,MiddleToTop,midHgating,delayGate,[timestamp(PostTime(1)) timestamp(PostTime(end))],divisionTime);
                p_peaks(1:size(p_temp,1),1:size(p_temp,2)) = double(p_temp);
                t_peaks(1:size(t_temp,1),1:size(t_temp,2)) = double(t_temp);
                combinedpeaks(1,1,:,:) = p_peaks;
                combinedpeaks(1,2,:,:) = t_peaks;
                h5write(fullfile(ndpathname,H5filename), peak_name, combinedpeaks, [double(scell) 1 1 1], [1 2 4 200]);

                clear p_params p_names t_params t_names param_names
                % calculating parameters from the detected Peaks
                [p_params,p_names] = peaksParam(p_temp);
                [t_params,t_names] = peaksParam(t_temp);
                
                
                % putting together parameter matrix
                param_mat = [param_mat;...
                    double(row),double(col),double(field),double(scell),double(phenotype),...
                    double(p_params),...
                    double(t_params),...
                    double(firstPeakDuration(t_peaks,timestamp(PostTime(1)),timestamp(PostTime(end)))),...
                    ];
                
                % putting together parameter name 
                
                param_names{1} = 'Row';
                param_names{2} = 'Column';
                param_names{3} = 'Field';
                param_names{4} = 'Cell Number';
                param_names{5} = 'Phenotype';
                oldsize = length(param_names);
                for i=1:length(p_names)
                    param_names{oldsize+i} = ['Pre-' p_names{i}];
                end
                oldsize = length(param_names);
                for i=1:length(t_names)
                    param_names{oldsize+i} = ['Post-' t_names{i}];
                end
                
                param_names{length(param_names)+1} = '1st peak delay';
                
                
            end
        end
        
        % if param_mat was assigned, write the new param_mat to the HDF5 file
        if numel(param_mat)>0
            fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
            
            if ~H5L.exists(fid,param_name,'H5P_DEFAULT')
                H5F.close(fid);
                display(['Initializing ' H5filename ':' param_name]);
            else
                %if the param_mat already exists, first delete the old param_mat
                H5L.delete(fid,param_name,'H5P_DEFAULT');
                display(['Overwriting ' H5filename ':' param_name]);
                H5F.close(fid);
            end
            
            % store parameters
            h5create(fullfile(ndpathname,H5filename), param_name, [size(param_mat,1), size(param_mat,2)], 'Datatype', 'double', 'ChunkSize', [size(param_mat,1), size(param_mat,2)], 'Deflate', 9);
            h5write(fullfile(ndpathname,H5filename), param_name, param_mat, [1 1], [size(param_mat,1) size(param_mat,2)]);
            
            for i = 1:size(param_mat,2)
                h5writeatt(fullfile(ndpathname,H5filename),param_name,['param' num2str(i)],param_names{i});
            end
            
        end
        display(['Finished calculating parameters of ' H5filename]);
        
    else
        display([signal_name ', ' timestamp_name ' or ' selectedcells_name ' does not exist']);
    end
else
    display([H5filename ' does not exist']);
end


function T = firstPeakDuration(peaks,refTime,endTime)
if ~isempty(peaks)
    firstPeakInd = find(peaks(4,:)==1,1,'first');
    if ~isempty(firstPeakInd)
        T = peaks(1,firstPeakInd) - refTime;
    else
        T = endTime;
    end
else
    T = endTime;
end


function [params, names]= peaksParam(peaks)

names{1} = 'Peak count';
names{2} = 'Mean(peak height)';
names{3} = 'SD(peak height)';
names{4} = 'Mean(peak duration)';
names{5} = 'SD(peak duration)';
names{6} = 'High state time';
names{7} = 'dH/dt, slope';
names{8} = 'dH/dt, R-square';
names{9} = 'dD/dt, slope';
names{10} = 'dD/dt, R-square';
names{11} = 'Mean(peak interval)';
names{12} = 'SD(peak interval)';


if ~isempty(peaks)

    params(1) = length(peaks(1,:));   % peak Count
    params(2) = mean(peaks(2,:));   % mean of peak height
    params(3) = std(peaks(2,:));      % sd of peak height
    params(4) = mean(peaks(3,:));   % mean of peak duration
    params(5) = std(peaks(3,:));      % sd of peak duration
    params(6) = sum(peaks(3,:));      % total time in high state
    params(7:8) = changerate(peaks(1,:),peaks(2,:)); % Height change rate: slope, rsquare
    params(9:10) = changerate(peaks(1,:),peaks(3,:)); % Duration change rate: slope, rsquare
    
    divisionPeakInd = find(peaks(4,:)==0);
    
    switch length(divisionPeakInd)
        case 1
            peakIntervals = [diff(peaks(1,1:divisionPeakInd(1)-1)) ...
                             diff(peaks(1,divisionPeakInd(1)+1:end))];
        case 2
            peakIntervals = [diff(peaks(1,1:divisionPeakInd(1)-1)) ...
                             diff(peaks(1,divisionPeakInd(1)+1:divisionPeakInd(2)-1)) ...
                             diff(peaks(1,divisionPeakInd(2)+1:end))];
        case 3
            peakIntervals = [diff(peaks(1,1:divisionPeakInd(1)-1)) ...
                             diff(peaks(1,divisionPeakInd(1)+1:divisionPeakInd(2)-1)) ...
                             diff(peaks(1,divisionPeakInd(2)+1:divisionPeakInd(3)-1)) ...
                             diff(peaks(1,divisionPeakInd(3)+1:end))];
        otherwise
            peakIntervals = diff(peaks(1,:));
    end
    if ~isempty(peakIntervals)
        params(11) = mean(peakIntervals); 
        params(12) = std(peakIntervals); 
    else
        params(11) = 0; 
        params(12) = 0; 
    end
else
    params(1:12) = [0 0 0 0 0 0 0 0 0 0 0 0];
end

function param = changerate(x,y)
if length(x)>=2
    [curve,gof] = fit( x', y', 'poly1');
    coeff = coeffvalues(curve);
    slope = coeff(1);
    if ~isnan(slope) 
        param(1,1) = slope;
    else
        param(1,1) = 0;
    end
    
    if ~isnan(gof.rsquare)
        param(1,2) = gof.rsquare;
    else
        param(1,2) = 0;
    end
    
else
    param = [0 0];
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

function output = findTruePeaks(xs,d_y,plotLog,ToTop,midHgating,delayGate,timewindows,divisionTime)
% output(1,:) = truePeak';
% output(2,:) = PeakHeight';
% output(3,:) = PeakDuration';
% output(4,:) = peakSelection;

[pks,locs] = findpeaks(d_y); %
options = statset('Display','off');
if isempty(pks) || length(locs)< 3
    output = [];
    return;
end

idx = kmeans(pks',3,'distance','sqEuclidean','emptyaction','drop','options',options,'Replicates',5);



matOrder(1,:) = [mean(pks(idx == 1)) mean(pks(idx == 2)) mean(pks(idx == 3))];
matOrder(2,:) = [std(pks(idx == 1)) std(pks(idx == 2)) std(pks(idx == 3))];
matOrder(3,:) = [1 2 3];
matOrder = sortrows(matOrder',1)';

TopP = locs(idx == matOrder(3,3));
BotP = locs(idx == matOrder(3,1));

if ToTop
    TopP = sort([TopP locs(idx == matOrder(3,2))]);
end

[pks2,locs2] = findpeaks(-d_y);
all_data = double([]);

all_data(1,:) = ([ones(1,length(pks)) -1*ones(1,length(pks2))]);
all_data(2,:) = ([pks -pks2]);
all_data(3,:) = ([locs locs2]);

sorted_data = sortrows(all_data',3)';
sorted_data(4,:) = diff([d_y(1) sorted_data(2,:)]);

idxn = kmeans(abs(sorted_data(4,:)),3,'distance','sqEuclidean','emptyaction','drop','options',options,'Replicates',5);

matOrder(1,:) = [mean(abs(sorted_data(4,idxn==1))) mean(abs(sorted_data(4,idxn==2))) mean(abs(sorted_data(4,idxn==3)))];
matOrder(2,:) = [1 2 3];
matOrder = sortrows(matOrder',1)';

TopD = find(idxn == matOrder(2,3));
BotD = find(idxn == matOrder(2,1));
median_BotD = median(abs(sorted_data(4,BotD)));
if ToTop
    midD = find(idxn == matOrder(2,2));
    
    for i = midD'
        if abs(sorted_data(4,i)) > midHgating% midHgating*median_BotD
            TopD = sort([TopD;i]);
        end
    end
end

PosTopDTail = [];
PosTopDTailH = [];
countind=1;

PosTopD = intersect(TopD, find(sorted_data(4,:)>0.04 & sorted_data(4,:)<0.4));
%median_PosTopD = median(sorted_data(4,PosTopD));
if ~isempty(PosTopD)
    for i=PosTopD'
        if i ~= length(sorted_data(4,:))
            decayH = 0;
            for j=(i+1):length(sorted_data(4,:))
                decayH = decayH+sorted_data(4,j);
                if -decayH > delayGate*sorted_data(4,i)
                    PosTopDTail(countind,1)  = sorted_data(3,j);
                    PosTopDTailH(countind,1) = sorted_data(4,i);
                    countind=countind+1;
                    break;
                end
            end
            if  -decayH < delayGate*sorted_data(4,i)
                PosTopDTail(countind,1)  = length(d_y); %sorted_data(3,end)+1;%
                PosTopDTailH(countind,1) = sorted_data(4,i);
                countind=countind+1;
            end
            
        else
            PosTopDTail(countind,1)  = length(d_y); %sorted_data(3,end)+1;%
            PosTopDTailH(countind,1) = sorted_data(4,i);
            countind=countind+1;
        end
    end
else
    output = [];
    return;
end

PosTopD_t = PosTopD(  xs(sorted_data(3,PosTopD))  >=timewindows(1)  &  xs(sorted_data(3,PosTopD))  <timewindows(2)  );
PosTopDTail = PosTopDTail(  xs(sorted_data(3,PosTopD))  >=timewindows(1)  &  xs(sorted_data(3,PosTopD))  <timewindows(2)  );
PosTopDTailH = PosTopDTailH(  xs(sorted_data(3,PosTopD))  >=timewindows(1)  &  xs(sorted_data(3,PosTopD))  <timewindows(2)  );
PosTopD = PosTopD_t;
countind = 1;
truePeak = [];
PeakHeight = [];
PeakDuration = [];
i=1;
while i<=length(PosTopDTail)
    
    temp = find(PosTopDTail == PosTopDTail(i));
    
    
    replicatedInd = temp(find( temp > i ));
    if isempty(replicatedInd)
        PeakStart = PosTopD(i);
        PeakTail = PosTopDTail(i);
        PHeight = PosTopDTailH(i);
        i=i+1;
    else
        PeakStart = PosTopD(i);
        PeakTail = PosTopDTail(i);
        PHeight = PosTopDTailH(i);
        i=max(replicatedInd)+1;
    end
    
    
    
    
    truePeak(countind)     = xs(sorted_data(3,PeakStart));
    PeakHeight(countind)   = PHeight;
    PeakDuration(countind) = xs(PeakTail) - truePeak(countind);
    countind=countind+1;

end
if isempty(truePeak)
    output = [];
    return;
end

% make sure that all peaks stay in the requested temporal windows
selected_truePeak = truePeak(truePeak+PeakDuration<=timewindows(2) &  truePeak>timewindows(1));
selected_PeakHeight = PeakHeight(truePeak+PeakDuration<=timewindows(2) &  truePeak>timewindows(1));
selected_PeakDuration = PeakDuration(truePeak+PeakDuration<=timewindows(2) &  truePeak>timewindows(1));
peakSelection = ones(1,length(selected_truePeak));
% get rid of peaks that stay in division time range
if ~isempty(divisionTime)
    for i = 1:length(divisionTime)
        frontPeakOverlapDivTime = find(selected_truePeak > divisionTime(i)-60 & selected_truePeak < divisionTime(i)+30 & selected_PeakDuration < 300);
        tailPeakOverlapDivTime = find(selected_truePeak+selected_PeakDuration > divisionTime(i)-60 & selected_truePeak+selected_PeakDuration < divisionTime(i)+30 & selected_PeakDuration < 300);
        DivisionWithinPeak = find(selected_truePeak<divisionTime(i)-60 & selected_truePeak+selected_PeakDuration > divisionTime(i)+30 & selected_PeakDuration < 300);
        divPeaks = union(union(frontPeakOverlapDivTime,tailPeakOverlapDivTime),DivisionWithinPeak);
        
        if ~isempty(divPeaks)
            nondivPeaks = setdiff(1:length(selected_truePeak),divPeaks);
            peakSelection(divPeaks) = 0;
            
        end
    end
end
if isempty(selected_truePeak)
    output = [];
    return;
else
    output(1,:) = selected_truePeak';
    output(2,:) = selected_PeakHeight';
    output(3,:) = selected_PeakDuration';
    output(4,:) = peakSelection;
end

if plotLog
    figure(2);
    
    subplot(2,1,2);
    plot(xs,d_y,'b-'); hold on;
    xlim([xs(1) xs(end)]);

    h = subplot(2,1,2);
    YLim = get(h,'YLim');
    
    for i=find(peakSelection==1)
        if selected_truePeak(i)+PeakDuration(i)<=timewindows(2) &&  selected_truePeak(i)>timewindows(1)
            rectangle('Position',[selected_truePeak(i),YLim(1),PeakDuration(i),YLim(2)-YLim(1)],...
                'FaceColor',[0.8 0.95 0.95],'EdgeColor','none','EraseMode','normal')
        else
            rectangle('Position',[selected_truePeak(i),YLim(1),PeakDuration(i),YLim(2)-YLim(1)],...
                'FaceColor',[0.9 0.9 0.9],'EdgeColor','none','EraseMode','normal')
        end
    end
    plot(xs,d_y,'b-'); hold on;
    plot(xs(TopP),d_y(TopP),'rv',xs(BotP),d_y(BotP),'gv');
    plot(xs(locs2),d_y(locs2),'k.');
    myylim = get(h,'Ylim');
    plot(timewindows',[myylim(2);myylim(2)],'ks-','MarkerFaceColor','k','Linewidth',5,'MarkerSize',8,'MarkerEdgeColor','none');
    if ~isempty(divisionTime)
        switch length(divisionTime)
            case 1
                plot([divisionTime(1)-60;divisionTime(1)+30],[myylim(1);myylim(1)],'r-','Linewidth',3);
            case 2
                plot([divisionTime(1)-60;divisionTime(1)+30],[myylim(1);myylim(1)],'r-','Linewidth',3);
                plot([divisionTime(2)-60;divisionTime(2)+30],[myylim(1);myylim(1)],'g-','Linewidth',3);
            case 3
                plot([divisionTime(1)-60;divisionTime(1)+30],[myylim(1);myylim(1)],'r-','Linewidth',3);
                plot([divisionTime(2)-60;divisionTime(2)+30],[myylim(1);myylim(1)],'g-','Linewidth',3);
                plot([divisionTime(3)-60;divisionTime(3)+30],[myylim(1);myylim(1)],'b-','Linewidth',3);
        end
    end
    hold off;
    
    subplot(2,1,1);
    stem(xs(sorted_data(3,PosTopD)),sorted_data(4,PosTopD),'r','filled'); hold on;
    stem(xs(PosTopDTail),PosTopDTailH','k','filled'); 
    xlim([xs(1) xs(end)]);
    hold off;
    %param_mat = [param_mat;double(scell) pvec1];
    clear all_data;
    
    pause(1);
    drawnow;
end

