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

%ndpathname = '/home/ss240/files/ImStor/sorger/data/NIC/Pat/02-03-2013/';
cellFateColor = [];
fateColor = [0 1 0;0.3 0.3 0.3;0.6 0.6 0.6;1 0 0;0 0 1];
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

            alldata = [alldata;param_mat(:,[13:25]+1)];
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
                cellFate =  [cellFate;param_mat(c_cell,27)];
                cellFateColor =  [cellFateColor;fateColor(param_mat(c_cell,27)+1,:)];
            end
            clear param_mat;
            
            
        else
            H5F.close(fid);
        end
        
end


w = 1./var(alldata);
[wcoeff,score,latent,tsquared,explained] = pca(alldata,'VariableWeights',w);

mapped_data = tsne(alldata,cellFate,score,30)
%mapped_data = compute_mapping(alldata,'SNE',3)
cGroup = 1:7;
testInd = [];
mycolor = jet(length(cGroup));
newColor = [];
cInds= [];
for i=1:length(cGroup)
    cInds = find(groupNo==cGroup(i));
    testInd = [testInd;cInds];
    for j=1:length(cInds)
        newColor=[newColor;mycolor(i,:)];
    end
end
testInd = sort(testInd);
figure(1)
scatter3(mapped_data(testInd,1),mapped_data(testInd,2),mapped_data(testInd,3),10,newColor,'o','fill'); hold on;

cGroup = setdiff(1:length(searchInd(1,:)),cGroup);
testInd = [];
mycolor = [0.3 0.3 0.3];
newColor = [];
cInds = [];
for i=1:length(cGroup)
    cInds = find(groupNo==cGroup(i));
    testInd = [testInd;cInds];
end
testInd = sort(testInd);
scatter3(mapped_data(testInd,1),mapped_data(testInd,2),mapped_data(testInd,3),5,[0.3 0.3 0.3],'fill'); hold off;


cGroup = 4:7;
testInd = [];
mycolor = jet(length(cGroup));
newColor = [];
for i=1:length(cGroup)
    cInds = find(groupNo==cGroup(i));
    testInd = [testInd;cInds];
    for j=1:length(cInds)
        newColor=[newColor;mycolor(i,:)];
    end
end
testInd = sort(testInd);
figure(2),scatter3(mapped_data(testInd,1),mapped_data(testInd,2),mapped_data(testInd,3),10,newColor,'o','fill'); hold on;

cGroup = setdiff(1:length(searchInd(1,:)),cGroup);
testInd = [];
mycolor = [0.3 0.3 0.3];
newColor = [];
cInds = [];
for i=1:length(cGroup)
    cInds = find(groupNo==cGroup(i));
    testInd = [testInd;cInds];
end
testInd = sort(testInd);
scatter3(mapped_data(testInd,1),mapped_data(testInd,2),mapped_data(testInd,3),5,[0.3 0.3 0.3],'fill'); hold off;
figure(3);
scatter3(mapped_data(:,1),mapped_data(:,2),mapped_data(:,3),10,cellFateColor,'fill'); hold off;
