function createClusterParams()
clc;
ndfilename = '02032013-r1.nd';
sourcefolder = 'c:\computation\02-03-2013';
%sourcefolder = '~/files/ImStor/sorger/data/NIC/Pat/02-03-2013/'
prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(sourcefolder,ndfilename);
outputsignalNo = 1;
sequenceNo = 2;
sites = [1:7 15:18];

if matlabpool('size') == 0
  matlabpool open;
end

parfor s=1:length(sites)
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
midHgating = 0.09; % x the median of lowest peak cluster
delayGate = 0.8; % fraction of height must decay to consider as peak tail

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
        c_peaks = [];
        
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
                ys=interp1(c_time,c_signal,xs);
                
                signal_mean=interp1(timestamp,signal_mean_t,xs);
                ys = ys-signal_mean;
                outTS = getTimeSeriesTrend(ys,'trendType',1);
                ys = outTS.dTS;
                p_peaks{scell} = findTruePeaks(xs,ys,0,MiddleToTop,midHgating,delayGate,[timestamp(PreTime(1)) timestamp(PreTime(end))],divisionTime);
                if isempty(p_peaks{scell})
                    p_peaks{scell} = [0;0;0];
                end
                
                %c_peaks{scell} = findTruePeaks(xs(PostTime(1):end),d_y(PostTime(1):end),showPlots,MiddleToTop,midHgating,delayGate);
                c_peaks{scell} = findTruePeaks(xs,ys,showPlots,MiddleToTop,midHgating,delayGate,[timestamp(PostTime(1)) timestamp(PostTime(end))],divisionTime);
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
            pause(1);
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

function output = findTruePeaks(xs,d_y,plotLog,ToTop,midHgating,delayGate,timewindows,divisionTime)
% output(1,:) = truePeak';
% output(2,:) = PeakHeight';
% output(3,:) = PeakDuration';

[pks,locs] = findpeaks(d_y); %
options = statset('Display','off');
if isempty(pks) || length(locs)< 3F
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
selected_truePeak = truePeak(truePeak+PeakDuration>=timewindows(1) &  truePeak<timewindows(2));
selected_PeakHeight = PeakHeight(truePeak+PeakDuration>=timewindows(1) &  truePeak<timewindows(2));
selected_PeakDuration = PeakDuration(truePeak+PeakDuration>=timewindows(1) &  truePeak<timewindows(2));

% get rid of peaks that stay in division time range

if ~isempty(divisionTime)
    for i = 1:length(divisionTime)
        frontPeakOverlapDivTime = find(selected_truePeak > divisionTime(i)-50 & selected_truePeak < divisionTime(i)+10 & selected_PeakDuration < 200);
        tailPeakOverlapDivTime = find(selected_truePeak+selected_PeakDuration > divisionTime(i)-50 & selected_truePeak+selected_PeakDuration < divisionTime(i)+10 & selected_PeakDuration < 200);
        DivisionWithinPeak = find(selected_truePeak<divisionTime(i)-50 & selected_truePeak+selected_PeakDuration > divisionTime(i)+10 & selected_PeakDuration < 200);
        divPeaks = union(union(frontPeakOverlapDivTime,tailPeakOverlapDivTime),DivisionWithinPeak);
        if ~isempty(divPeaks)
            nondivPeaks = setdiff(1:length(selected_truePeak),divPeaks);
            selected_truePeak = selected_truePeak(nondivPeaks);
            selected_PeakHeight = selected_PeakHeight(nondivPeaks);
            selected_PeakDuration = selected_PeakDuration(nondivPeaks);
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
end

if plotLog
    figure(2);
    
    subplot(2,1,2);
    plot(xs,d_y,'b-'); hold on;
    xlim([xs(1) xs(end)]);

    h = subplot(2,1,2);
    YLim = get(h,'YLim');
    
    for i=1:length(selected_truePeak)
        if selected_truePeak(i)+PeakDuration(i)>=timewindows(1) &&  selected_truePeak(i)<timewindows(2)
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
                plot([divisionTime(1)-50;divisionTime(1)+10],[myylim(1);myylim(1)],'r-','Linewidth',3);
            case 2
                plot([divisionTime(1)-50;divisionTime(1)+10],[myylim(1);myylim(1)],'r-','Linewidth',3);
                plot([divisionTime(2)-50;divisionTime(2)+10],[myylim(1);myylim(1)],'g-','Linewidth',3);
            case 3
                plot([divisionTime(1)-50;divisionTime(1)+10],[myylim(1);myylim(1)],'r-','Linewidth',3);
                plot([divisionTime(2)-50;divisionTime(2)+10],[myylim(1);myylim(1)],'g-','Linewidth',3);
                plot([divisionTime(3)-50;divisionTime(3)+10],[myylim(1);myylim(1)],'b-','Linewidth',3);
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
    
    
    drawnow;
end



%function detrended_y = detrendingSignal(xs,yori,plotLog,ToTop,midHgating,delayGate,signal_mean)


% 
% [pks,locs] = findpeaks(ys); %
% options = statset('Display','off');
% if isempty(pks) || length(locs)< 3
%     output = [];
%     return;
% end
% 
% idx = kmeans(pks',3,'distance','sqEuclidean','emptyaction','drop','options',options,'Replicates',5);
% 
% 
% 
% matOrder(1,:) = [mean(pks(idx == 1)) mean(pks(idx == 2)) mean(pks(idx == 3))];
% matOrder(2,:) = [std(pks(idx == 1)) std(pks(idx == 2)) std(pks(idx == 3))];
% matOrder(3,:) = [1 2 3];
% matOrder = sortrows(matOrder',1)';
% 
% TopP = locs(idx == matOrder(3,3));
% BotP = locs(idx == matOrder(3,1));
% 
% if ToTop
%     TopP = sort([TopP locs(idx == matOrder(3,2))]);
% end
% 
% [pks2,locs2] = findpeaks(-ys);
% all_data = double([]);
% 
% all_data(1,:) = ([ones(1,length(pks)) -1*ones(1,length(pks2))]);
% all_data(2,:) = ([pks -pks2]);
% all_data(3,:) = ([locs locs2]);
% 
% sorted_data = sortrows(all_data',3)';
% sorted_data(4,:) = diff([ys(1) sorted_data(2,:)]);
% 
% idxn = kmeans(abs(sorted_data(4,:)),3,'distance','sqEuclidean','emptyaction','drop','options',options,'Replicates',5);
% 
% matOrder(1,:) = [mean(abs(sorted_data(4,idxn==1))) mean(abs(sorted_data(4,idxn==2))) mean(abs(sorted_data(4,idxn==3)))];
% matOrder(2,:) = [1 2 3];
% matOrder = sortrows(matOrder',1)';
% 
% TopD = find(idxn == matOrder(2,3));
% BotD = find(idxn == matOrder(2,1));
% median_BotD = median(abs(sorted_data(4,BotD)));
% if ToTop
%     midD = find(idxn == matOrder(2,2));
%     
%     for i = midD'
%         if abs(sorted_data(4,i)) > midHgating% midHgating*median_BotD
%             TopD = sort([TopD;i]);
%         end
%     end
% end
% 
% PosTopDTail = [];
% PosTopDTailH = [];
% countind=1;
% 
% PosTopD = intersect(TopD, find(sorted_data(4,:)>0.04 & sorted_data(4,:)<0.4));
% %median_PosTopD = median(sorted_data(4,PosTopD));
% if ~isempty(PosTopD)
% for i=PosTopD'
%     if i ~= length(sorted_data(4,:))
%         decayH = 0;
%         for j=(i+1):length(sorted_data(4,:))
%             decayH = decayH+sorted_data(4,j);
%             if -decayH > delayGate*sorted_data(4,i)
%                 PosTopDTail(countind,1)  = sorted_data(3,j);
%                 PosTopDTailH(countind,1) = sorted_data(4,i);
%                 countind=countind+1;
%                 break;
%             end
%         end
%         if  -decayH < delayGate*sorted_data(4,i)
%             PosTopDTail(countind,1)  = length(ys); %sorted_data(3,end)+1;%
%             PosTopDTailH(countind,1) = sorted_data(4,i);
%             countind=countind+1;
%         end
% 
%     else
%         PosTopDTail(countind,1)  = length(ys); %sorted_data(3,end)+1;%
%         PosTopDTailH(countind,1) = sorted_data(4,i);
%         countind=countind+1;
%     end
% end
% else
%     output = [];
%     return;
% end
% 
% 
% SigP = [];
% tailP = [];
% i=1;
% while i<=length(PosTopDTail)
%     
%     temp = find(PosTopDTail == PosTopDTail(i));
%     
%     
%     replicatedInd = temp(find( temp > i ));
%     if isempty(replicatedInd)
%         PeakStart = PosTopD(i);
%         PeakTail = PosTopDTail(i);
%         i=i+1;
%     else
%         PeakStart = PosTopD(i);
%         PeakTail = PosTopDTail(i);
%         i=max(replicatedInd)+1;
%     end
%     
%     SigP = [SigP sorted_data(3,PeakStart)];
%     tailP = [tailP PeakTail];
% end
% 
% 
% ind=1;
% for i=1:length(BotP)
%     FloorP(ind) = BotP(i);
%     ind = ind+1;
% end
% 
% FloorP = setdiff(union(setdiff(FloorP,SigP),tailP),length(xs));
% NOISE = [xs(FloorP);ys(FloorP)];
% SIGNAL = [xs(SigP);ys(SigP)];
% 
% pp = csaps(NOISE(1,:),NOISE(2,:));
% y_sp = fnval(pp,xs) ; 
% x_sp_2half = xs(round(length(xs)/2):FloorP(end));
% y_sp_2half = y_sp(round(length(xs)/2):FloorP(end));
% final_y = interp1(x_sp_2half,y_sp_2half,xs(end),'linear','extrap');
% 
% new_x = [xs(1) NOISE(1,:) xs(end)];
% new_y = [ys(1) NOISE(2,:) final_y];
% w = 10*ones(size(new_x));
% w(new_y(1:(length(FloorP))) > (mean(new_y)+1*std(new_y))) = 1e-4;
% 
% pp2 = csaps(new_x,new_y,0,[],w);
% y_sp2 = fnval(pp2,xs) ; 
% detrended_y = ys-y_sp2;
% 

% figure(1);
% plot(xs,ys,'-k'); hold on;
% %plot(xs(end),final_y,'ok');
% %plot(SIGNAL(1,:),SIGNAL(2,:),'og',NOISE(1,:),NOISE(2,:),'xr'); 
% %fnplt(pp2,'r--');
% plot(xs,detrended_y,'-g');
% hold off;
% drawnow;
% 
% 

% [pks,locs] = findpeaks(ys);
% options = statset('Display','off');
% 
% %obj = gmdistribution.fit(pks',2,'Options',options,'Replicates',3,'SharedCov',false);
% %idx = cluster(obj,pks');
% if length(pks)<50
%     detrended_y = ys;
%     return;
% end
% idx = kmeans(pks',2,'distance','sqEuclidean','emptyaction','drop','options',options,'Replicates',5);
% 
% if abs(median(pks(idx == 1)) - median(pks(idx == 2))) < 0.09
%     idx = ones(size(idx));
% end
% temp_idx = idx;
% if median(-pks(idx == 1)) < median(-pks(idx == 2))
%     temp_idx = ones(size(idx));
%     temp_idx(idx == 1) = 2;
% end
% idx = temp_idx;
% 
% [pksn,locsn] = findpeaks(-ys);
% 
% %obj2 = gmdistribution.fit(-(pksn'),2,'Options',options,'Replicates',3);
% %idxn = cluster(obj2,-(pksn'));
% idxn = kmeans(-(pksn'),2,'distance','sqEuclidean','emptyaction','drop','options',options,'start','uniform','Replicates',5);
% temp_idxn = idxn;
% if median(-pksn(idxn == 1)) > median(-pksn(idxn == 2))
%     temp_idxn = ones(size(idxn));
%     temp_idxn(idxn == 1) = 2;
% end
% idxn = temp_idxn;
% 
% 
% if abs(median(pksn(idxn == 1)) - median(pksn(idxn == 2))) < 0.09
%     idxn = ones(size(idxn));
% end
% 
% 
% all_data(1,:) = ([ones(1,length(pks)) -1*ones(1,length(pksn))]);
% all_data(2,:) = ([idx' idxn']);
% all_data(3,:) = ([pks -pksn]);
% all_data(4,:) = ([locs locsn]);
% sorted_data = sortrows(all_data',4)';
% sorted_data(5,:) = diff([ys(1) sorted_data(3,:)]); %true peak height
% 
% peakDist = abs(sorted_data(5,:)) ;
% %obj3 = gmdistribution.fit(peakDist',2,'Options',options);
% %idx3 = cluster(obj3,peakDist');
% idx3 = kmeans(peakDist',2,'distance','sqEuclidean','emptyaction','drop','options',options,'Replicates',5);
% 
% temp_idx3 = idx3;
% if median(peakDist(idx3==1)) > median(peakDist(idx3==2))
%     temp_idx3 = ones(size(idx3));
%     temp_idx3(idx3 == 1) = 2;
% end
% idx3 = temp_idx3;
% SmallH = find(idx3 == 1)   ;
% smallPosH = intersect(SmallH , find(sorted_data(5,:)>0) )  ;
% largeH = find(idx3 == 2);
% LargeNegH = setdiff(intersect(largeH , find(sorted_data(5,:)<0) ) , find(sorted_data(2,:)==2) )  ;
% LargePosH = intersect(largeH , find(sorted_data(5,:)>0) )  ;
% selectedLargeNegH = LargeNegH (  abs(sorted_data(5,LargeNegH)) > HeightThres  ) ;
% unselectedlargeH = setdiff(largeH,selectedLargeNegH);
% foundList = [];
% SelectedSmallH = [];
% segmentNo = 0;
% for i = 1:length(selectedLargeNegH)
%     if ~isempty(find(selectedLargeNegH(i)+1 == SmallH,1))
%         foundList = [foundList selectedLargeNegH(i)];
%         for j=(selectedLargeNegH(i)+1):length(idx3)
%             if  ~ isempty(find(j == SmallH,1) )
%                 SelectedSmallH = [SelectedSmallH j];
%                 segmentNo=segmentNo+1;
%             else
%                 break;
%             end
%         end
%     end
% end
% 
% if smallPosHLoc
%     
%     for i = 1:length(smallPosH)
%         if ~isempty(find(smallPosH(i)+1 == SmallH,1))
%             for j=(smallPosH(i)+1):length(idx3)
%                 if  ~ isempty(find(j == SmallH,1) )
%                     SelectedSmallH = [SelectedSmallH j];
%                     segmentNo=segmentNo+1;
%                 else
%                     break;
%                 end
%             end
%         end
%     end
%     
% end
% 
% truesmallH = find(idx3 == 1);
% 
% if truesmallH(1) == 1
%     for j=1:length(sorted_data(4,:))
%         if  ~ isempty(find(j == truesmallH,1) )
%             SelectedSmallH = [SelectedSmallH j];
%             segmentNo=segmentNo+1;
%         else
%             break;
%         end
%     end
% end
% 
% SelectedSmallH = sort(SelectedSmallH);
% 
% baseInd = union(SelectedSmallH,selectedLargeNegH);
% if size(baseInd,2) > 1
%     baseInd = baseInd';
% end
% LowWeightInd = union(setdiff(find(idx3 == 1),baseInd),unselectedlargeH);
% if size(LowWeightInd,2)>1
%     LowWeightInd = LowWeightInd';
% end
% spList=[];
% spList(:,1) = [baseInd;LowWeightInd];
% spList(:,2) = [100*ones(size(baseInd)); 0.01*ones(size(LowWeightInd))];
% spList = sortrows(spList,1);
% 
% spX = xs(sorted_data(4,spList(:,1)));
% spY = sorted_data(3,spList(:,1));
% spW = spList(:,2)';
% 
% if ~isempty(LargePosH)
% if LargePosH(1) == 1
%     spX = [xs(1),spX];
%     spY = [ys(1),spY];
%     spW = [1  ,spW];
%     
% else
%     spX = [xs(1),spX];
%     spY = [ys(1),spY];
%     spW = [0.01  ,spW];
% end
% end
% 
% 
% if length(sorted_data(3,:)) == smallPosH(end)
%     for j=length(sorted_data(3,:)):-1:1
%         if isempty(find(j == truesmallH,1))
%             break;
%         end
%     end
%     
%     substractingH = sorted_data(5,j);
%     if substractingH>0
%         spX = [spX xs(end)];
%         spY = [spY ys(end)- mode(abs(sorted_data(5,largeH)))];
%         spW = [spW 30];
%     else
%         spX = [spX xs(end)];
%         spY = [spY ys(end)];
%         spW = [spW 0.01];
%     end
%     
% else
%     
%     spX = [spX xs(end)];
%     spY = [spY ys(end)];
%     spW = [spW 0.01];
% end
% %in_knots = abs(diff(xs(sorted_data(4,baseInd)))) > 100;
% %[xs(1) xs(1) xs(sorted_data(4,baseInd(in_knots)))  xs(end) xs(end) ]
% sp_base = spap2(noknots,2,spX,spY,spW);
% y_sp = fnval(sp_base,xs) ;
% 
% detrended_y = ys-y_sp;
% 
% if plotLog
%     figure(1);
%     subplot(4,1,1);
%     plot(xs,ys);hold on;
%     
%     plot(xs(locs(idx == 1)),pks(idx == 1),'rv',xs(locs(idx == 2)),pks(idx == 2),'gv',xs(locs(idx == 3)),pks(idx == 3),'bv')
%     plot(xs(locsn(idxn == 1)),-pksn(idxn == 1),'r^',xs(locsn(idxn == 2)),-pksn(idxn == 2),'g^',xs(locsn(idxn == 3)),-pksn(idxn == 3),'b^')
%     hold off;
%     
%     subplot(4,1,2);
%     stem(xs(sorted_data(4,sorted_data(2,:)==1)),sorted_data(5,sorted_data(2,:)==1),'filled','r'); hold on;
%     stem(xs(sorted_data(4,sorted_data(2,:)==2)),sorted_data(5,sorted_data(2,:)==2),'filled','g'); hold off;
%     
%     subplot(4,1,3);
%     stem(xs(sorted_data(4,idx3==1)),peakDist(idx3==1),'r'); hold on;
%     stem(xs(sorted_data(4,SmallH)),peakDist(SmallH),'r','filled');
%     stem(xs(sorted_data(4,idx3==2)),peakDist(idx3==2),'g');
%     stem(xs(sorted_data(4,selectedLargeNegH)),peakDist(selectedLargeNegH),'g','filled');
%     plot(xs(sorted_data(4,foundList)),peakDist(foundList),'yx');
%     plot(xs(sorted_data(4,LargePosH)),peakDist(LargePosH),'bx');
%     plot(xs(sorted_data(4,SelectedSmallH)),peakDist(SelectedSmallH),'ko');
%     hold off;
%     
%     subplot(4,1,4);
%     plot(xs,ys);hold on;
%     plot(xs(sorted_data(4,baseInd)),sorted_data(3,baseInd),'.r');
%     plot(xs(sorted_data(4,LowWeightInd)),sorted_data(3,LowWeightInd),'.g');
%     plot(xs,y_sp,'k-');
%     
%     hold off;
%     
% end


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

