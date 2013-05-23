function createClusterParams()
clc;
ndfilename = '02032013-r1.nd';
sourcefolder = 'c:\computation\02-03-2013';
%sourcefolder = '~/files/ImStor/sorger/data/NIC/Pat/02-03-2013/'
prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(sourcefolder,ndfilename);
outputsignalNo = 3;
sequenceNo = 2;
sites = [15:21];

if matlabpool('size') == 0
  matlabpool open;
end

parfor site=sites
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
HeightThres = 0.01;
mypeakThres = 0.001;
smallPosHLoc= 1;
MiddleToTop = 1;
showPlots = 1;
noknots = 1;
midHgating = 6; % x the median of lowest peak cluster
delayGate = 0.3; % fraction of height must decay to consider as peak tail

PreTime = 1:25;
PostTime = 32:225;

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
        cellpath = cell(size(cellpath_mat,3),1);
        sisterList = cell(size(cellpath_mat,3),1);
        
        for tp=1:size(cellpath_mat,3)
            cellpath{tp} = cellpath_mat(:,:,tp);
            sisterList{tp} = sisterList_mat(:,:,tp);
        end
        
        param_mat = [];
        c_peaks = [];
        for scell = setdiff(selected_cells',[]); %%%%
            
            
            
            %if find(scell==selected_cells,1)
            PosTime = find(signal(:,scell)~=0);
            c_signal = signal(PosTime,scell);
            c_time   = timestamp(PosTime);
            if cellpath{end}(scell,1) > 0
                if sisterList{end}(scell,1) ~= -1
                    phenotype = 2; %dividing
                else
                    phenotype = 1; %non-dividing
                end
            elseif cellpath{end}(scell,1) == -2
                phenotype = 0;     %died
                display('cell dead');
            end
               
            
            if numel(PosTime)>minimumSignalSize
                display([H5filename 'cell:' num2str(scell)]);
                xs = c_time(1):7:timestamp(PosTime(end));
                ys=interp1(c_time,c_signal,xs);
                
                if median(ys) > 0.9
                    d_y = detrendingSignal(xs,ys,showPlots,HeightThres,mypeakThres,smallPosHLoc,noknots);
                else
                    outTS = getTimeSeriesTrend(ys,'trendType',1);
                    d_y = outTS.dTS;
                end
                % output(1,:) - truePeak
                % output(2,:) - PeakHeight
                % output(3,:) - PeakDuration
                p_peaks{scell} = findTruePeaks(xs(PreTime),d_y(PreTime),0,MiddleToTop,midHgating,delayGate);
                if isempty(p_peaks{scell})
                    p_peaks{scell} = [0;0;0];
                end
                %outTS = getTimeSeriesTrend(d_y(PostTime),'trendType',1);

                %d_y2 = detrendingSignal(xs(PostTime),d_y(PostTime),0,HeightThres,mypeakThres,smallPosHLoc,noknots);
                c_peaks{scell} = findTruePeaks(xs(PostTime(1):end),d_y(PostTime(1):end),showPlots,MiddleToTop,midHgating,delayGate);
                if isempty(c_peaks{scell})
                    c_peaks{scell} = [0;0;0];
                end
                
                param_mat = [param_mat;...
                    double(scell),...
                    median(ys(PreTime)),...            % 1 median of pre-signal
                    median(ys(PostTime(1):end)),...           % 2 median of post-signal
                    length(p_peaks{scell}(1,:)),...    % 3 peak Count
                    median(p_peaks{scell}(2,:)),...    % 4 median peak height
                    median(p_peaks{scell}(3,:)),...    % 5 median peak duration
                    sum(p_peaks{scell}(3,:)),...       % 6 total time in high state
                    changerate(p_peaks{scell}(1,:),p_peaks{scell}(2,:)),... % 7 8 Height change rate: slope rsquare
                    changerate(p_peaks{scell}(1,:),p_peaks{scell}(3,:)),... % 9 10 Duration change rate: slope rsquare
                    length(c_peaks{scell}(1,:)),...    % 11 peak Count
                    median(c_peaks{scell}(2,:)),...    % 12 median peak height
                    median(c_peaks{scell}(3,:)),...    % 13 median peak duration
                    sum(c_peaks{scell}(3,:)),...       % 14 total time in high state
                    changerate(c_peaks{scell}(1,:),c_peaks{scell}(2,:)),... % 15 16 Height change rate: slope rsquare
                    changerate(c_peaks{scell}(1,:),c_peaks{scell}(3,:)),... % 17 18 Duration change rate: slope rsquare
                    phenotype,...
                    ];


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
        end
        display(['Finished calculating parameters of ' H5filename]);
        
    else
        display([signal_name ', ' timestamp_name ' or ' selectedcells_name ' does not exist']);
    end
else
    display([H5filename ' does not exist']);
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

function output = findTruePeaks(xs,d_y,plotLog,ToTop,midHgating,delayGate)
% output(1,:) = truePeak';
% output(2,:) = PeakHeight';
% output(3,:) = PeakDuration';

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
        if abs(sorted_data(4,i)) > midHgating*median_BotD
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
    output = []
    return;
end
output(1,:) = truePeak';
output(2,:) = PeakHeight';
output(3,:) = PeakDuration';


if plotLog
    figure(2);
    subplot(2,1,1);
    stem(xs(sorted_data(3,PosTopD)),sorted_data(4,PosTopD),'r','filled'); hold on;
    %stem(xs(sorted_data(3,BotD)),sorted_data(4,BotD),'b'); hold on;
    stem(xs(PosTopDTail),PosTopDTailH,'k','filled'); hold on;
    xlim([xs(1) xs(end)]);
    hold off;
    %param_mat = [param_mat;double(scell) pvec1];
    clear all_data;
    
    subplot(2,1,2);
    plot(xs,d_y,'b-'); hold on;
    xlim([xs(1) xs(end)]);

    h = subplot(2,1,2);
    YLim = get(h,'YLim');
    
    for i=1:length(truePeak)
        rectangle('Position',[truePeak(i),YLim(1),PeakDuration(i),YLim(2)-YLim(1)],...
            'FaceColor',[0.9 0.9 0.9],'EdgeColor','none','EraseMode','normal')
    end
    plot(xs,d_y,'b-'); hold on;
    plot(xs(TopP),d_y(TopP),'rv',xs(BotP),d_y(BotP),'gv');
    plot(xs(locs2),d_y(locs2),'k.')
    hold off;
    drawnow;
end

function detrended_y = detrendingSignal(xs,ys,plotLog,HeightThres,mypeakThres,smallPosHLoc,noknots)
[pks,locs] = findpeaks(ys);
options = statset('Display','off');

%obj = gmdistribution.fit(pks',2,'Options',options,'Replicates',3,'SharedCov',false);
%idx = cluster(obj,pks');
if length(pks)<50
    detrended_y = ys;
    return;
end
idx = kmeans(pks',2,'distance','sqEuclidean','emptyaction','drop','options',options,'Replicates',5);

if abs(median(pks(idx == 1)) - median(pks(idx == 2))) < 0.09
    idx = ones(size(idx));
end
temp_idx = idx;
if median(-pks(idx == 1)) < median(-pks(idx == 2))
    temp_idx = ones(size(idx));
    temp_idx(idx == 1) = 2;
end
idx = temp_idx;

[pksn,locsn] = findpeaks(-ys);

%obj2 = gmdistribution.fit(-(pksn'),2,'Options',options,'Replicates',3);
%idxn = cluster(obj2,-(pksn'));
idxn = kmeans(-(pksn'),2,'distance','sqEuclidean','emptyaction','drop','options',options,'start','uniform','Replicates',5);
temp_idxn = idxn;
if median(-pksn(idxn == 1)) > median(-pksn(idxn == 2))
    temp_idxn = ones(size(idxn));
    temp_idxn(idxn == 1) = 2;
end
idxn = temp_idxn;


if abs(median(pksn(idxn == 1)) - median(pksn(idxn == 2))) < 0.09
    idxn = ones(size(idxn));
end


all_data(1,:) = ([ones(1,length(pks)) -1*ones(1,length(pksn))]);
all_data(2,:) = ([idx' idxn']);
all_data(3,:) = ([pks -pksn]);
all_data(4,:) = ([locs locsn]);
sorted_data = sortrows(all_data',4)';
sorted_data(5,:) = diff([ys(1) sorted_data(3,:)]);

peakDist = abs(sorted_data(5,:)) ;
%obj3 = gmdistribution.fit(peakDist',2,'Options',options);
%idx3 = cluster(obj3,peakDist');
idx3 = kmeans(peakDist',2,'distance','sqEuclidean','emptyaction','drop','options',options,'Replicates',5);

temp_idx3 = idx3;
if median(peakDist(idx3==1)) > median(peakDist(idx3==2))
    temp_idx3 = ones(size(idx3));
    temp_idx3(idx3 == 1) = 2;
end
idx3 = temp_idx3;
SmallH = find(idx3 == 1)   ;
smallPosH = intersect(SmallH , find(sorted_data(5,:)>0) )  ;
largeH = find(idx3 == 2);
LargeNegH = setdiff(intersect(largeH , find(sorted_data(5,:)<0) ) , find(sorted_data(2,:)==2) )  ;
LargePosH = intersect(largeH , find(sorted_data(5,:)>0) )  ;
selectedLargeNegH = LargeNegH (  abs(sorted_data(5,LargeNegH)) > HeightThres  ) ;
unselectedlargeH = setdiff(largeH,selectedLargeNegH);
foundList = [];
SelectedSmallH = [];
segmentNo = 0;
for i = 1:length(selectedLargeNegH)
    if ~isempty(find(selectedLargeNegH(i)+1 == SmallH,1))
        foundList = [foundList selectedLargeNegH(i)];
        for j=(selectedLargeNegH(i)+1):length(idx3)
            if  ~ isempty(find(j == SmallH,1) )
                SelectedSmallH = [SelectedSmallH j];
                segmentNo=segmentNo+1;
            else
                break;
            end
        end
    end
end

if smallPosHLoc
    
    for i = 1:length(smallPosH)
        if ~isempty(find(smallPosH(i)+1 == SmallH,1))
            for j=(smallPosH(i)+1):length(idx3)
                if  ~ isempty(find(j == SmallH,1) )
                    SelectedSmallH = [SelectedSmallH j];
                    segmentNo=segmentNo+1;
                else
                    break;
                end
            end
        end
    end
    
end

truesmallH = find(idx3 == 1);

if truesmallH(1) == 1
    for j=1:length(sorted_data(4,:))
        if  ~ isempty(find(j == truesmallH,1) )
            SelectedSmallH = [SelectedSmallH j];
            segmentNo=segmentNo+1;
        else
            break;
        end
    end
end

SelectedSmallH = sort(SelectedSmallH);

baseInd = union(SelectedSmallH,selectedLargeNegH);
if size(baseInd,2) > 1
    baseInd = baseInd';
end
LowWeightInd = union(setdiff(find(idx3 == 1),baseInd),unselectedlargeH);
if size(LowWeightInd,2)>1
    LowWeightInd = LowWeightInd';
end
spList=[];
spList(:,1) = [baseInd;LowWeightInd];
spList(:,2) = [100*ones(size(baseInd)); 0.01*ones(size(LowWeightInd))];
spList = sortrows(spList,1);

spX = xs(sorted_data(4,spList(:,1)));
spY = sorted_data(3,spList(:,1));
spW = spList(:,2)';

if ~isempty(LargePosH)
if LargePosH(1) == 1
    spX = [xs(1),spX];
    spY = [ys(1),spY];
    spW = [1  ,spW];
    
else
    spX = [xs(1),spX];
    spY = [ys(1),spY];
    spW = [0.01  ,spW];
end
end


if length(sorted_data(3,:)) == smallPosH(end)
    for j=length(sorted_data(3,:)):-1:1
        if isempty(find(j == truesmallH,1))
            break;
        end
    end
    
    substractingH = sorted_data(5,j);
    if substractingH>0
        spX = [spX xs(end)];
        spY = [spY ys(end)- mode(abs(sorted_data(5,largeH)))];
        spW = [spW 30];
    else
        spX = [spX xs(end)];
        spY = [spY ys(end)];
        spW = [spW 0.01];
    end
    
else
    
    spX = [spX xs(end)];
    spY = [spY ys(end)];
    spW = [spW 0.01];
end
%in_knots = abs(diff(xs(sorted_data(4,baseInd)))) > 100;
%[xs(1) xs(1) xs(sorted_data(4,baseInd(in_knots)))  xs(end) xs(end) ]
sp_base = spap2(noknots,2,spX,spY,spW);
y_sp = fnval(sp_base,xs) ;

detrended_y = ys-y_sp;

if plotLog
    figure(1);
    subplot(4,1,1);
    plot(xs,ys);hold on;
    
    plot(xs(locs(idx == 1)),pks(idx == 1),'rv',xs(locs(idx == 2)),pks(idx == 2),'gv',xs(locs(idx == 3)),pks(idx == 3),'bv')
    plot(xs(locsn(idxn == 1)),-pksn(idxn == 1),'r^',xs(locsn(idxn == 2)),-pksn(idxn == 2),'g^',xs(locsn(idxn == 3)),-pksn(idxn == 3),'b^')
    hold off;
    
    subplot(4,1,2);
    stem(xs(sorted_data(4,sorted_data(2,:)==1)),sorted_data(5,sorted_data(2,:)==1),'filled','r'); hold on;
    stem(xs(sorted_data(4,sorted_data(2,:)==2)),sorted_data(5,sorted_data(2,:)==2),'filled','g'); hold off;
    
    subplot(4,1,3);
    stem(xs(sorted_data(4,idx3==1)),peakDist(idx3==1),'r'); hold on;
    stem(xs(sorted_data(4,SmallH)),peakDist(SmallH),'r','filled');
    stem(xs(sorted_data(4,idx3==2)),peakDist(idx3==2),'g');
    stem(xs(sorted_data(4,selectedLargeNegH)),peakDist(selectedLargeNegH),'g','filled');
    plot(xs(sorted_data(4,foundList)),peakDist(foundList),'yx');
    plot(xs(sorted_data(4,LargePosH)),peakDist(LargePosH),'bx');
    plot(xs(sorted_data(4,SelectedSmallH)),peakDist(SelectedSmallH),'ko');
    hold off;
    
    subplot(4,1,4);
    plot(xs,ys);hold on;
    plot(xs(sorted_data(4,baseInd)),sorted_data(3,baseInd),'.r');
    plot(xs(sorted_data(4,LowWeightInd)),sorted_data(3,LowWeightInd),'.g');
    plot(xs,y_sp,'k-');
    
    hold off;
    
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

