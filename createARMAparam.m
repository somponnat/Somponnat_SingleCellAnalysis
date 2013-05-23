function createARMAparam()
clc;
ndfilename = '02032013-r1.nd';
sourcefolder = 'c:\computation\02-03-2013';
%sourcefolder = '~/files/ImStor/sorger/data/NIC/Pat/02-03-2013/'
prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(sourcefolder,ndfilename);
outputsignalNo = 3;
sequenceNo = 2;
sites = [4];

if matlabpool('size') == 0
   matlabpool open;
end

for site=sites
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

minimumSignalSize = 170;

H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
signal_name = ['/field' num2str(field)  '/outputsignal' num2str(outputsignalNo)];
timestamp_name = ['/field' num2str(field) '/timestamp' num2str(outputsignalNo)];
selectedcells_name = ['/field' num2str(field) '/selectedcells'];
param_name = ['/field' num2str(field)  '/clusterparams' num2str(outputsignalNo)];

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
        param_mat = [];

        for scell = selected_cells'
            %if find(scell==selected_cells,1)
                PosTime = find(signal(40:end,scell));
                c_signal = signal(PosTime,scell);
                c_time   = timestamp(PosTime);
                
                if numel(c_time)>minimumSignalSize && numel(c_signal)>minimumSignalSize
                    display([H5filename 'cell:' num2str(scell)]);
                    ys=interp1(c_time,c_signal,c_time(1):7:c_time(end));
                    outTS = getTimeSeriesTrend(ys,'trendType',4);
                    
                    D1 = LagOp({1,-1},'Lags',[0,1]);
                    
                    xData = 100:length(outTS.dTS);
                    yData =  filter(D1,outTS.dTS(100:end))';
                    %Y=outTS.dTS(100:end)';
                    model = arima(4,0,1);
                    options = optimset('fmincon');
                    options = optimset(options,'Display','off','Diagnostics','off');
                    [fit,VarCov,LogL,info] = estimate(model,yData,'options',options,'print',false);
                    [coefficients, ~] = toCellArray(LagOp(fit.AR)/LagOp(fit.MA));
                    pvec1 = cepstralCoff(20,cell2mat(coefficients));
                    figure(1);
                    plot(yData,'b'); hold on;
                    plot(simulate(fit,length(yData)),'r'); hold off;
                    title(num2str(pvec1));
                    
                    %                 % Set up fittype and options.
                    %                 ft = fittype( 'poly8' );
                    %                 opts = fitoptions( ft );
                    %                 opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf];
                    %                 opts.Upper = [Inf Inf Inf Inf Inf Inf Inf Inf Inf];
                    
                    %
                    %
                    %                 mydata2 = iddata(outTS.trend',[],7,'Tstart',0);
                    %                 assignin('base','ys',outTS.trend')
                    %
                    %                 Opt = arxOptions;
                    %                 arx4 = arx(mydata,[4], Opt);
                    %                 pvec2  = getpvec(arx4);
                    
                    param_mat = [param_mat;double(scell) pvec1];
                end
            %end
            
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
        end
        display(['Finished calculating parameters of ' H5filename]);
        
    else
        display([signal_name ', ' timestamp_name ' or ' selectedcells_name ' does not exist']);
    end
else
    display([H5filename ' does not exist']);
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

