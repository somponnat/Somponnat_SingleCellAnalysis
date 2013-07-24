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
selected_truePeak = truePeak(truePeak+PeakDuration>=timewindows(1) &  truePeak<timewindows(2));
selected_PeakHeight = PeakHeight(truePeak+PeakDuration>=timewindows(1) &  truePeak<timewindows(2));
selected_PeakDuration = PeakDuration(truePeak+PeakDuration>=timewindows(1) &  truePeak<timewindows(2));
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
    
    
    drawnow;
end
