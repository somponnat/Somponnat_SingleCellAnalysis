function createClusterParams()
clc;
ndfilename = '02032013-r1.nd';
sourcefolder = 'C:\computation\02-03-2013';
%sourcefolder = '~/files/ImStor/sorger/data/NIC/Pat/02-03-2013/'
prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(sourcefolder,ndfilename);
outputsignalNo = 1;
sequenceNo = 2;
sites = [29];

if matlabpool('size') == 0
  matlabpool open;
end

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



function cal_clusterparam(row,col,field,ndpathname,outputsignalNo,sequenceNo)
minimumSignalSize = 100;
MiddleToTop = 1;
showPlots = 0;
midHgating = 0.1; % x the median of lowest peak cluster
delayGate = 0.7; % fraction of height must decay to consider as peak tail

PreTime = 1:25;
PostTime = 40:244;

H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
signal_name = ['/field' num2str(field)  '/outputsignal' num2str(outputsignalNo)];
timestamp_name = ['/field' num2str(field) '/timestamp' num2str(outputsignalNo)];
selectedcells_name = ['/field' num2str(field) '/selectedcells'];
param_name = ['/field' num2str(field)  '/clusterparams' num2str(outputsignalNo)];
cellpath_name = ['/field' num2str(field) '/cellpath'];
sisterList_name = ['/field' num2str(field) '/sisterList'];

warning off;

if exist(fullfile(ndpathname,H5filename),'file')
    fileattrib(fullfile(ndpathname,H5filename),'+w');
    fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if H5L.exists(fid,signal_name,'H5P_DEFAULT') && ...
            H5L.exists(fid,timestamp_name,'H5P_DEFAULT') && ...
            H5L.exists(fid,selectedcells_name,'H5P_DEFAULT')
        H5F.close(fid);
        signalinfo = h5info(fullfile(ndpathname,H5filename), signal_name);
        startind = double([1 1 sequenceNo]);
        countind = [signalinfo.Dataspace.Size(1) signalinfo.Dataspace.Size(2) 1];
        signal = permute(h5read(fullfile(ndpathname,H5filename),signal_name,startind, countind),[2 1 3]);
        timestamp = h5read(fullfile(ndpathname,H5filename),timestamp_name);
        selected_cells = h5read(fullfile(ndpathname,H5filename),selectedcells_name);
        %outputsignal_name = h5readatt(fullfile(ndpathname,H5filename),signal_name,['signal' num2str(sequenceNo)]);
        cellpathinfo = h5info(fullfile(ndpathname,H5filename), cellpath_name);
        sisterListinfo = h5info(fullfile(ndpathname,H5filename), sisterList_name);

        cellpath_mat = h5read(fullfile(ndpathname,H5filename),cellpath_name,[1 1 1], [cellpathinfo.Dataspace.Size(1) cellpathinfo.Dataspace.Size(2) cellpathinfo.Dataspace.Size(3)]);
        sisterList_mat = h5read(fullfile(ndpathname,H5filename),sisterList_name,[1 1 1], [sisterListinfo.Dataspace.Size(1) sisterListinfo.Dataspace.Size(2) sisterListinfo.Dataspace.Size(3)]);
        
        param_mat = [];
        
        signal_mean_t = smooth(mean(signal,2));
        
        for scell = setdiff(selected_cells',[]); %%%%
            
            if cellpath_mat(scell,1,end) > 0
                if     sisterList_mat(scell,1,end) ~= -1 && sisterList_mat(scell,2,end) == -1 && sisterList_mat(scell,3,end) == -1
                    phenotype = 2; %dividing once
                elseif sisterList_mat(scell,1,end) ~= -1 && sisterList_mat(scell,2,end) ~= -1 && sisterList_mat(scell,3,end) == -1
                    phenotype = 3; %dividing twice
                elseif sisterList_mat(scell,1,end) ~= -1 && sisterList_mat(scell,2,end) ~= -1 && sisterList_mat(scell,3,end) ~= -1
                    phenotype = 4; %dividing three times
                else
                    phenotype = 1; %non-dividing
                end
            elseif cellpath_mat(scell,1,end) == -2
                phenotype = 0;     %died
                display('cell dead');
            end
            
            %find division time
            divisionTime = [];
            switch phenotype
                case 2
                    divisionTime(1) = timestamp(find(sisterList_mat(scell,1,:) ~= -1,1,'first'));
                case 3
                    divisionTime(1) = timestamp(find(sisterList_mat(scell,1,:) ~= -1,1,'first'));
                    divisionTime(2) = timestamp(find(sisterList_mat(scell,2,:) ~= -1,1,'first'));
                case 4
                    divisionTime(1) = timestamp(find(sisterList_mat(scell,1,:) ~= -1,1,'first'));
                    divisionTime(2) = timestamp(find(sisterList_mat(scell,2,:) ~= -1,1,'first'));
                    divisionTime(3) = timestamp(find(sisterList_mat(scell,3,:) ~= -1,1,'first'));
                case {0,1}
                    divisionTime = [];
            end
            
            PosTime = find(signal(:,scell)~=0);
            c_signal = signal(PosTime,scell);
            c_time   = timestamp(PosTime);            
               
            if numel(PosTime)>minimumSignalSize
                display([H5filename 'cell:' num2str(scell)]);
                xs = c_time(1):7:timestamp(PosTime(end));
                ys_ori=interp1(c_time,c_signal,xs);
                ys=ys_ori;
                
                signal_mean=interp1(timestamp,signal_mean_t,xs);
                ys = ys_ori-signal_mean;
                outTS = getTimeSeriesTrend(ys,'trendType',1);
                ys = outTS.dTS;
                p_peaks{scell} = findTruePeaks(xs,ys,0,MiddleToTop,midHgating,delayGate,[timestamp(PreTime(1)) timestamp(PreTime(end))],divisionTime);
                
                t_peaks{scell} = findTruePeaks(xs,ys,showPlots,MiddleToTop,midHgating,delayGate,[timestamp(PostTime(1)) timestamp(PostTime(end))],divisionTime);
                
                clear p_params p_names t_params t_names param_names
                [p_params,p_names] = peaksParam(p_peaks{scell});
                [t_params,t_names] = peaksParam(t_peaks{scell});
                
                % putting together parameter matrix
                param_mat = [param_mat;...
                    double(row),double(col),double(field),double(scell),...
                    double(p_params),...
                    double(t_params),...
                    double(firstPeakDuration(t_peaks{scell},timestamp(PostTime(1)))),...
                    double(phenotype),...
                    ];
                
                % putting together parameter name struct
                param_names{1} = 'Row';
                param_names{2} = 'Column';
                param_names{3} = 'Field';
                param_names{4} = 'Cell Number';
                oldsize = length(param_names);
                for i=1:length(p_names)
                    param_names{oldsize+i} = ['Pre-' p_names{i}];
                end
                oldsize = length(param_names);
                for i=1:length(t_names)
                    param_names{oldsize+i} = ['Post-' t_names{i}];
                end
                
                param_names{length(param_names)+1} = '1st peak delay';
                param_names{length(param_names)+1} = 'Phenotype';
            end
        end
        if numel(param_mat)>0
            fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
            
            if ~H5L.exists(fid,param_name,'H5P_DEFAULT')
                H5F.close(fid);
                display(['Initializing ' H5filename ':' param_name]);
            else
                H5L.delete(fid,param_name,'H5P_DEFAULT');
                display(['Overwriting ' H5filename ':' param_name]);
                H5F.close(fid);
            end
            
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

function T = firstPeakDuration(peaks,refTime)
if ~isempty(peaks)
    firstPeakInd = find(peaks(4,:)==1,1,'first');
    if ~isempty(firstPeakInd)
        T = peaks(1,firstPeakInd) - refTime;
    else
        T = -1;
    end
else
    T = -1;
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

