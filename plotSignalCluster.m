close all; clear all; clc;
alldata = [];
names = [];
groupColor = [];
groupNo = [];
cellFate = [];
field = 1;
outputsignalNo =1;
searchInd(1,:) = [1 1 1 1 1 1 1 3 3 3 3 ];
searchInd(2,:) = [3 4 5 6 7 8 9 3 4 5 6 ];

%searchInd(1,:) = [3 3 3 3 ];
%searchInd(2,:) = [3 4 5 6 ];

%searchInd(1,:) = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4  5  6  7  8];
%searchInd(2,:) = [3 5 7 8 9 3 5 7 8 9 3 4 5 6 7 3 4 5 6 7 10 10 10 10];

mycolor = [hsv(5);hsv(5);hsv(5);hsv(5);hsv(5)];%size(searchInd,2));
mycolor2 = hsv(3);
ndpathname = 'c:\computation\02-03-2013';
noCluster = 20;
%ndpathname = '/home/ss240/files/ImStor/sorger/data/NIC/Pat/02-03-2013/';

data_index = 1;
for i = 1:size(searchInd,2)
        row = searchInd(1,i);
        col = searchInd(2,i);
        
        H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
        param_name = ['/field' num2str(field)  '/clusterparams' num2str(outputsignalNo)];
        signal_name = ['/field' num2str(field)  '/outputsignal' num2str(outputsignalNo)];
        timestamp_name = ['/field' num2str(field) '/timestamp' num2str(outputsignalNo)];
        sisterList_name = ['/field' num2str(field) '/sisterList'];
        
        fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
        if H5L.exists(fid,param_name,'H5P_DEFAULT')
            H5F.close(fid);
            paraminfo = h5info(fullfile(ndpathname,H5filename), param_name);
            startind = double([1 1]);
            countind = [paraminfo.Dataspace.Size(1) paraminfo.Dataspace.Size(2)];
            param_mat = double(h5read(fullfile(ndpathname,H5filename),param_name,startind, countind));

            alldata = [alldata;param_mat(:,[11:18]+1)];
            groupNo = [groupNo;i*ones(size(param_mat,1),1)];
            
            signalinfo = h5info(fullfile(ndpathname,H5filename), signal_name);
            %sisterListinfo = h5info(fullfile(handles.SourceF,H5filename), sisterList_name);
            %sisterList = h5read(fullfile(handles.SourceF,H5filename),sisterList_name,[1 1 1],...
            %                 [sisterListinfo.Dataspace.Size(1) sisterListinfo.Dataspace.Size(2) sisterListinfo.Dataspace.Size(3)]);
            startind = double([1 1 2]);
            countind = [signalinfo.Dataspace.Size(1) signalinfo.Dataspace.Size(2) 1];
            signal = permute(h5read(fullfile(ndpathname,H5filename),signal_name,startind, countind),[2 1 3]);
            timestamp = h5read(fullfile(ndpathname,H5filename),timestamp_name);
            for c_cell=1:size(param_mat,1)
                PosTime = find(signal(:,param_mat(c_cell,1)));
                c_signal = signal(PosTime,param_mat(c_cell,1));
                c_time   = timestamp(PosTime);
                ys=interp1(c_time,c_signal,0:7:1568,'spline');
                
                originData(data_index,:) = ys;
                data_index=data_index+1;
                names{length(names)+1} = ['r' num2str(row) 'c' num2str(col) '_' num2str(param_mat(c_cell,1))];
                groupColor = [groupColor;mycolor(i,:)];
                cellFate =  [cellFate;param_mat(c_cell,20)];
            end
            clear param_mat;
            
            
        else
            H5F.close(fid);
        end
        
end


w = 1./var(alldata);
[wcoeff,score,latent,tsquared,explained] = pca(alldata,'VariableWeights',w);

%c3 = wcoeff(:,1:3)
coefforth = diag(sqrt(w))*wcoeff;

cscores = zscore(alldata)*coefforth;

figure(1)

gscatter(score(:,1),score(:,2),cellFate,'bgrkm','ssxoo',3,'on','1st Principal Component','2nd Principal Component');

xlabel('1st Principal Component');
ylabel('2nd Principal Component');

figure(14)
interestedG = 10;
gscatter(score(groupNo==interestedG,1),score(groupNo==interestedG,2),cellFate(groupNo==interestedG),mycolor2,'x',3,'on','1st Principal Component','2nd Principal Component');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');

plotG = 6;
figure(15);set(gcf,'Position',[500 500 400 400]);
myX = 1;
myY = 2;
gscatter(score(groupNo~=plotG,myX),score(groupNo~=plotG,myY),groupNo(groupNo~=plotG),[0.6 0.6 0.6],'o',2.5,'off');
hold on;
gscatter(score(groupNo==plotG-1,myX),score(groupNo==plotG-1,myY),groupNo(groupNo==plotG-1),[0 1 0],'o',2.5,'off');
gscatter(score(groupNo==plotG-2,myX),score(groupNo==plotG-2,myY),groupNo(groupNo==plotG-2),[0 0 0],'o',2.5,'off');
gscatter(score(groupNo==plotG+1,myX),score(groupNo==plotG+1,myY),groupNo(groupNo==plotG+1),[0 0 1],'o',2.5,'off');
gscatter(score(groupNo==plotG,myX),score(groupNo==plotG,myY),groupNo(groupNo==plotG),[1 0 0],'o',2.5,'off','1st Principal Component','2nd Principal Component');
axis equal;

figure(16);
mycolor3=jet(5);
for i=11:15
    plotG=i;plot(0:7:1568,mean(originData(groupNo==plotG,:),1),'color',mycolor3(i-10,:)); hold on
end


%gname(names)

figure(2);
pareto(explained);
xlabel('Principal Component');
ylabel('Variance Explained (%)');

[st2,index] = sort(tsquared,'descend'); % sort in descending order
extreme = index(1:10);

figure(3);
biplot(coefforth(:,1:2),'scores',score(:,1:2),'varlabels',num2str((1:size(alldata,2))'));


figure(6);
for i=1:4%size(wcoeff,2)
    subplot(2,2,i),bar(coefforth(:,i));
    title(['Component' num2str(i)]);
    
end
figure(5);
biplot(coefforth(:,1:3),'scores',score(:,1:3),'varlabels',num2str((1:size(alldata,2))'),'Marker','x');

hold off;
%axis([-.26 0.8 -.51 .51 -.61 .81]);
%view([30 40]);
figure(4), 
scatter3(score(:,1),score(:,2),score(:,3),22,groupColor);
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
zlabel('3rd Principal Component');


Z = linkage(score(:,1:3),'ward','euclidean');
T = cluster(Z,'maxclust',noCluster);
%T = kmeans(score(:,1:2),noCluster,'replicates',4);
% 
% 



mycolor = hsv(noCluster);

figure(7)
scatter3(score(:,1),score(:,2),score(:,3),22,mycolor(T,:),'x');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
zlabel('3rd Principal Component');
figure(10)
gscatter(score(:,1),score(:,2),T,mycolor,'osv',3,'on','1st Principal Component','2nd Principal Component');

%figure(8);

%D = pdist(score(:,1:3),'euclidean');
%leafOrder = optimalleaforder(Z,D);
%[H,T2,outperm] = dendrogram(Z,0,'Labels',num2str(T),'Orientation','left','ColorThreshold',60); % ,'Reorder',leafOrder

%figure(9);
%load MyColormaps
%heatmap(mat2gray(originData(outperm,:),[0.9 1.5]),[],T,[],...
%   'ColorLevels',128,'Colormap',mycmap1,'FontSize',1);



binning = [];
for i=1:size(searchInd,2)
    binSize=[];
    for j=1:noCluster
        binning{i,j} = find(groupNo==i & T==j);
        binSize(j) = numel(binning{i,j});
    end
    figure(20); subplot(3,7,i);
    bar(1:noCluster,binSize);
    title(['Row:' num2str(searchInd(1,i)) ', Col:' num2str(searchInd(2,i))]);
    xlim([0 noCluster+1]);%ylim([0 0.8]);
end
clear refGroupList myList1 taken_myList1 n x temp sizeTaken newList

refGroup = 1;
refGroupList = T(groupNo==refGroup);
uniqueCluster = unique(refGroupList);
temp = zeros(length(uniqueCluster),2);
for i=1:length(uniqueCluster)
    temp(i,1) = uniqueCluster(i);
    temp(i,2) = length(find(refGroupList==uniqueCluster(i)));
end

sizeTaken =length(uniqueCluster);
myList1 = flipud(sortrows(temp,2));
if size(myList1,1)<sizeTaken
    taken_myList1 = myList1(1:size(myList1,1),1)';
else
    taken_myList1 = myList1(1:sizeTaken,1)';
end

clear refGroupList taken_myList2  temp sizeTaken2 

refGroup = 7;

refGroupList = T(groupNo==refGroup);
uniqueCluster = unique(refGroupList);
temp2 = zeros(length(uniqueCluster),2);
for i=1:length(uniqueCluster)
    temp2(i,1) = uniqueCluster(i);
    temp2(i,2) = length(find(refGroupList==uniqueCluster(i)));
end

sizeTaken2 =length(uniqueCluster);
myList2 = flipud(sortrows(temp2,2));


if size(myList1,1)<sizeTaken2
    taken_myList2 = myList2(1:size(myList2,1),1)';
else
    taken_myList2 = myList2(1:sizeTaken2,1)';
end
clear newList;
newList = zeros(1,noCluster);
FrontI = 1;
BackI = noCluster;
for i=1:max([length(taken_myList1) length(taken_myList2)])
    if i<=length(taken_myList1) && isempty(find(newList==taken_myList1(i),1))
        newList(FrontI) = taken_myList1(i);
        taken_myList2 = taken_myList2(taken_myList2~=taken_myList1(i));
        FrontI=FrontI+1;
    end
    
    if i<=length(taken_myList2) && isempty(find(newList==taken_myList2(i),1))
        newList(BackI) = taken_myList2(i);
        taken_myList1 = taken_myList1(taken_myList1~=taken_myList2(i));
        BackI=BackI-1;
    end
    
end
emptyInd = find(newList==0);
emptyValue = setdiff(1:noCluster,newList);
for i=1:length(emptyInd)
    newList(emptyInd(i)) = emptyValue(i);
end

k = 1;
for i=1:size(searchInd,2)
    
    countInd=1;
    binSize=[];
    for j=newList
        binSize(countInd) = numel(binning{i,j});
        countInd=countInd+1;
    end
    
    figure(20); subplot(3,7,k);
    k = k+1;
    h=bar(1:noCluster,binSize/sum(binSize));
    set(gca,'XTickLabel',[]);
    %text(1:noCluster,0.2*ones(1,noCluster),str2num(newList'));
    
    xlim([0 noCluster+1]); ylim([0 0.7]);
    title(['Row:' num2str(searchInd(1,i)) ', Col:' num2str(searchInd(2,i))]);
end

% [clusterInfo, pointToClusterMap] = MeanShiftClustering(score(:,1:2),0.1);
% plotMeanShiftResult(score(:,1:2),clusterInfo,pointToClusterMap);
% binning = [];
% 
% noC = size(clusterInfo,2)
% for i=1:size(searchInd,2)
%     binSize=[];
%     for j=1:noC
%         binning{i,j} = find(groupNo==i & pointToClusterMap==j);
%         binSize(j) = numel(binning{i,j});
%     end
%     figure(11); subplot(4,5,i);
%     bar(1:noC,binSize);
%     title(['Row:' num2str(searchInd(1,i)) ', Col:' num2str(searchInd(2,i))]);
%     %xlim([0 noC+1]);%ylim([0 0.8]);
% end
% clear myList1 taken_myList1 n x temp sizeTaken newList1 newList2
% 
% refGroup = 1;
% 
% [n x] = hist(pointToClusterMap(groupNo==refGroup),unique(pointToClusterMap(groupNo==refGroup))');
% temp(:,1) = x';
% temp(:,2) = n';
% 
% sizeTaken =round(length(unique(pointToClusterMap(groupNo==refGroup))));
% myList1 = flipud(sortrows(temp,2));
% if size(myList1,1)<sizeTaken
%     taken_myList1 = myList1(1:size(myList1,1),1)';
% else
%     taken_myList1 = myList1(1:sizeTaken,1)';
% end
% 
% clear myList2 taken_myList2 n2 x2 temp2 sizeTaken2 
% 
% refGroup = 5;
% 
% [n2 x2] = hist(pointToClusterMap(groupNo==refGroup),unique(pointToClusterMap(groupNo==refGroup))');
% temp2(:,1) = x2';
% temp2(:,2) = n2';
% 
% sizeTaken2 =round(length(unique(pointToClusterMap(groupNo==refGroup))));
% myList2 = flipud(sortrows(temp2,2));
% if size(myList1,1)<sizeTaken2
%     taken_myList2 = myList2(1:size(myList2,1),1)';
% else
%     taken_myList2 = myList2(1:sizeTaken2,1)';
% end
% newList1 = zeros(1,noC);
% newList2 = zeros(1,noC);
% FrontI = 1;
% BackI = 1;
% for i=1:max([length(taken_myList1) length(taken_myList2)])
%     
%     if i<=length(taken_myList1)
%         
%         if isempty(find(newList2==taken_myList1(i),1))
%             newList1(FrontI) = taken_myList1(i);
%             FrontI=FrontI+1;
% %         elseif (find(newList1==taken_myList1(i)) <= find(newList2==taken_myList1(i),1))  && temp(taken_myList1(i)==temp(:,1),2)/sum(n1) > temp2(taken_myList1(i)==temp2(:,1),2)/sum(n2)
% %             
% %             newList2(newList2==taken_myList1(i),1)=0;
% %             newList2 = newList2(newList2~=0);
% %             newList1(FrontI) = taken_myList1(i);
% %             FrontI=FrontI+1;
% %             BackI=BackI-1;
%          end
%         
%     end
%     
%     if i<=length(taken_myList2)
%         if isempty(find(newList1==taken_myList2(i),1)) 
%             newList2(BackI) = taken_myList2(i);
%             BackI=BackI+1;
% %         elseif (find(newList2==taken_myList2(i)) <= find(newList1==taken_myList2(i),1)) && temp2(taken_myList2(i)==temp2(:,1),2)/sum(n2) > temp(taken_myList2(i)==temp(:,1),2)/sum(n1)
% %             
% %             newList1(newList1==taken_myList2(i),1)=0;
% %             newList1 = newList1(newList1~=0);
% %             
% %             newList2(BackI) = taken_myList2(i);
% %             BackI=BackI+1;
% %            FrontI=FrontI-1;
%          end
%         
%     end
%     
% end
% newList2 = fliplr(newList2);
% newList = zeros(1,noC);
% newList(newList1~=0) = newList1(newList1~=0);
% newList(newList2~=0) = newList2(newList2~=0);
% 
% %emptyInd = sort(find(newList==0));
% % 
% % refGroup = 14;
% % clear taken_myList3 n3 x3 temp3  
% % [n3 x3] = hist(pointToClusterMap(groupNo==refGroup),unique(pointToClusterMap(groupNo==refGroup))');
% % collectL = [];
% % for i=1:length(x3)
% %     if ~isempty(find(x3(i) == intersect(x3,emptyInd),1))
% %         collectL = [collectL i];
% %     end
% % end
% % 
% % temp3(:,1) = x3(collectL)';
% % temp3(:,2) = n3(collectL)';
% % taken_myList3 = flipud(sortrows(temp3,2))';
% % 
% % newList(emptyInd(1):(emptyInd(1)+length(taken_myList3)-1)) = taken_myList3(1,:)
% 
% newList(newList==0) = setdiff(1:noC,find(newList~=0));
% 
% k = 1;
% for i=1:size(searchInd,2)
%     
%     countInd=1;
%     binSize=[];
%     for j=newList
%         binSize(countInd) = numel(binning{i,j});
%         countInd=countInd+1;
%     end
%     
%     figure(11); subplot(4,5,k);
%     k = k+1;
%     h=bar(1:noC,binSize/sum(binSize));
%     set(gca,'XTickLabel',[]);
%     %text(1:noCluster,0.2*ones(1,noCluster),str2num(newList'));
%     
%     xlim([0 noC+1]); ylim([0 0.2]);
%     title(['Row:' num2str(searchInd(1,i)) ', Col:' num2str(searchInd(2,i))]);
% end