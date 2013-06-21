function PhenotypeAnalysis()
SourceF = 'C:\computation\02-03-2013\';
%SourceF = '~/files/ImStor/soger/data/NIC/Pat/02-03-2013';
ndfilename = '02032013-r1.nd';

site = 4;


field = 1;
[notp stagePos stageName channelnames] = readndfile(fullfile(SourceF,ndfilename));
if notp==-1
    return;
end

first_tp = 1;
last_tp = notp;
prefix = ndfilename(1:(end-3));
%fileformat = [prefix '_%s_s' num2str(site) '_t%g.TIF']
tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
if ~isempty(tokens)
    row = tokens{1}{1};
    col = tokens{1}{2};
else
    row =1;
    col =1;
end



H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
cellpath_name = ['/field' num2str(field) '/cellpath'];
sisterList_name = ['/field' num2str(field) '/sisterList'];
bg_name = ['/field' num2str(field) '/bg'];
timestamp_name = ['/field' num2str(field) '/timestamp1'];

cellpathinfo = h5info(fullfile(SourceF,H5filename), cellpath_name);
sisterListinfo = h5info(fullfile(SourceF,H5filename), sisterList_name);
bginfo = h5info(fullfile(SourceF,H5filename), bg_name);
cellpath_mat = h5read(fullfile(SourceF,H5filename),cellpath_name,[1 1 1], [cellpathinfo.Dataspace.Size(1) cellpathinfo.Dataspace.Size(2) cellpathinfo.Dataspace.Size(3)]);
sisterList_mat = h5read(fullfile(SourceF,H5filename),sisterList_name,[1 1 1], [sisterListinfo.Dataspace.Size(1) sisterListinfo.Dataspace.Size(2) sisterListinfo.Dataspace.Size(3)]);
bg_mat = h5read(fullfile(SourceF,H5filename),bg_name,[1 1 1], [bginfo.Dataspace.Size(1) bginfo.Dataspace.Size(2) bginfo.Dataspace.Size(3)]);
timestamp = double(h5read(fullfile(SourceF,H5filename),timestamp_name));
cellpath = cell(last_tp,1);
sisterList = cell(last_tp,1);


for tp=first_tp:size(cellpath_mat,3)
    cellpath{tp} = cellpath_mat(:,:,tp);
    sisterList{tp} = sisterList_mat(:,:,tp);

end
NoSurv = [];
NoEmer = [];
NoDead = [];
Dead_p1 = [];
NoDiv = [];
wDiv = [];
withDiv_frac = [];
noDiv_frac = [];
assignin('base','cellpath',cellpath);
for tp=1:last_tp
    c_cellpath = cellpath{tp};
    c_sisterList = sisterList{tp};
    survPop = find(c_cellpath(:,1)>0);
    survPop1 = find(cellpath{1}(:,1)>0);
    wSiPop = find(c_sisterList(:,1)~=-1);    %find(sisterList(cellNo,1,:)~=-1 & sisterList(cellNo,2,:)~=-1 & sisterList(cellNo,3,:)==-1,1,'first');
    noSiPop = find(c_sisterList(:,1)==-1);
    emerPop = find(c_cellpath(:,1)==-1);     %find(cellpath_mat(cellNo,1,:)==-2 & cellpath_mat(cellNo,2,:)==-2,1,'first');
    deadPop = find(c_cellpath(:,1)==-2);
    NoSurv = [NoSurv length(survPop)];
    NoEmer = [NoEmer length(emerPop)];
    NoDead = [NoDead length(deadPop)];
    Dead_p1 = [Dead_p1 length(intersect(survPop1,deadPop))];
    NoDiv = [NoDiv length(intersect(survPop1,intersect(survPop,noSiPop)))];
    wDiv = [wDiv length(intersect(survPop1,intersect(survPop,wSiPop)))];
    
    withDiv_frac = [withDiv_frac length(intersect(survPop,wSiPop))+length(intersect(deadPop,wSiPop))];
    noDiv_frac =   [noDiv_frac   length(intersect(survPop,noSiPop))+length(intersect(deadPop,noSiPop))];
    
end
tp = last_tp;
%for tp=1:last_tp
figure(1);
plot(timestamp(1:tp)/60,100.*Dead_p1(1:tp)./NoSurv(1),'b');hold on;
plot(timestamp(1:tp)/60,100.*NoDiv(1:tp)./NoSurv(1),'g');
plot(timestamp(1:tp)/60,100.*wDiv(1:tp)./NoSurv(1),'r');
plot(timestamp(1:tp)/60,100.*NoSurv(1:tp)./NoSurv(1),'k'); hold off;
%    drawnow;
%end


text(timestamp(last_tp)/60,100.*NoSurv(last_tp)./NoSurv(1),['Net:' num2str(100.*NoSurv(last_tp)./NoSurv(1),'%6.1f') '%'],'VerticalAlignment','bottom','HorizontalAlignment','right');
text(timestamp(last_tp)/60,100.*NoDiv(last_tp)./NoSurv(1),['Surviving ori pop-quiescent:' num2str(100.*NoDiv(last_tp)./NoSurv(1),'%6.1f') '%'],'VerticalAlignment','bottom','HorizontalAlignment','right');
text(timestamp(last_tp)/60,100.*wDiv(last_tp)./NoSurv(1),['Surviving ori pop-divided:' num2str(100.*wDiv(last_tp)./NoSurv(1),'%6.1f') '%'],'VerticalAlignment','bottom','HorizontalAlignment','right');

text(timestamp(last_tp)/60,100.*Dead_p1(last_tp)./NoSurv(1),['Dead ori pop:' num2str(100.*Dead_p1(last_tp)./NoSurv(1),'%6.1f') '%'],'VerticalAlignment','top','HorizontalAlignment','right');
figure(1); ylabel('Change in cell number (%)'); xlabel('Time (hours)');
title(['R' num2str(row) 'C' num2str(col) ' Start:' num2str(NoSurv(first_tp)) ' End:' num2str(NoSurv(last_tp)) ' Non-dividing:' num2str(NoDiv(last_tp)) ' Dead:' num2str(NoDead(last_tp))]);


figure(2);
tol = 3e-4;
mytotal = log2((NoSurv+NoDead)./NoSurv(1));
dividing = mytotal.*withDiv_frac./(noDiv_frac+withDiv_frac);
quiescent = mytotal.*noDiv_frac./(noDiv_frac+withDiv_frac);
assignin('base','time',timestamp(1:last_tp)/60)
assignin('base','mytotal',mytotal)
assignin('base','dividing',dividing)
assignin('base','quiescent',quiescent)
smooth_total = spaps(timestamp(1:last_tp)/60,mytotal(1:last_tp),tol);
smooth_dividing = spaps(timestamp(1:last_tp)/60,dividing(1:last_tp),tol);
smooth_quiescent = spaps(timestamp(1:last_tp)/60,quiescent(1:last_tp),tol);
%fnplt(smooth_total,'g'); 

%fnplt(smooth_dividing,'b');
%fnplt(smooth_quiescent,'r');
plot(timestamp(1:last_tp)/60,mytotal,'g'); hold on;
plot(timestamp(1:last_tp)/60,dividing,'b');
plot(timestamp(1:last_tp)/60,quiescent,'r'); hold off;

text(timestamp(last_tp)/60,mytotal(last_tp),['Total'],'VerticalAlignment','bottom','HorizontalAlignment','right','Color','g');
text(timestamp(last_tp)/60,dividing(last_tp),['Dividing'],'VerticalAlignment','bottom','HorizontalAlignment','right','Color','b');
text(timestamp(last_tp)/60,quiescent(last_tp),['Quiescent'],'VerticalAlignment','bottom','HorizontalAlignment','right','Color','r');
ylabel('Population doublings'); xlabel('Time (hours)');
title(['R' num2str(row) 'C' num2str(col) ' Start:' num2str(NoSurv(first_tp)) ' End:' num2str(NoSurv(last_tp)) ' Non-dividing:' num2str(NoDiv(last_tp)) ' Dead:' num2str(NoDead(last_tp))]);

function [notp stagePos stageName waveName] = readndfile(filename)
% Search for number of string matches per line.
notp=-1;
stagePos = [];
stageName = [];
waveName = [];


if exist(filename,'file')
    fid = fopen(filename);
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
            stagename  = regexp(tline, '(?<="Stage\d+", ")\w+(?=")', 'match');
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