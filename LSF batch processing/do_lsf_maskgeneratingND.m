function do_lsf_maskgeneratingND()
clear all;
% Define information about input images and necessary parameters-----------
ndfilename = '01262013-r1.nd';
templateCH = 4;
cellsize = 40;
celldilatesize = 4;
%sourcefolder = 'c:\computation\01-26-2013\s22';
sourcefolder = '/files/ImStor/sorger/data/NIC/Pat/01-26-2013';
%------------------------------------------------
currentF = pwd;
cd(sourcefolder);
prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(ndfilename);
cd(currentF);
tps = [1 notp];
%sites = [1:4 7:8 13:20];
sites = [5:6 9:12 21:24];

jobmgr = findResource('scheduler', 'type', 'lsf');
jobmgr.ClusterMatlabRoot = '/opt/matlab';
jobmgr.SubmitArguments = '-q long -W 30:00 -R "rusage[matlab_dc_lic=1]"';
job = jobmgr.createJob();


for i=1:length(sites)
    site = sites(i);
    fileformat = [prefix '_%s_s' num2str(site) '_t%g.TIF'];
    %tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)f(?<field>\d+)','tokens');
    tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
    
    if ~isempty(tokens)
        row = str2num(tokens{1}{1});
        col = str2num(tokens{1}{2});
    else
        row = site;
        col = 1;
    end
    
    %field = tokens{1}{3};
    field = 1;
    plane = 1;
    job.createTask(@maskgen_commandline, 0, ...
        {3,sourcefolder,row,col,field,plane,templateCH,tps,fileformat,channelnames,cellsize,celldilatesize});
end

job.submit();

function maskgen_commandline(filetype,targetfolder,row,col,field,plane,templateCH,tps,fileformat,channelnames,cellsize,celldilatesize)

currentF = pwd;
cd(targetfolder);

% Load celltrack
celltrackOUTfilename = ['celltrackOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(templateCH)];
maskOUTfilename = ['maskOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(templateCH)];

if ~exist([celltrackOUTfilename '.mat'],'file') 
    display([celltrackOUTfilename '.mat does not exist.']);
    return
end
% file exist, load the file
load(celltrackOUTfilename,'cellpath','sisterList','bg');
if isempty(cellpath) || isempty(bg) || isempty(sisterList) %#ok<*NODEF>
    display([celltrackOUTfilename '.mat does not contain necessary variables.']);
    return
end

firsttp = tps(1);
lasttp = tps(end);

selected_cells = [];
save(maskOUTfilename,'selected_cells');
templateinfo = loadimageinfo(filetype,fileformat,[row col field plane templateCH],firsttp,channelnames);
imwidth = templateinfo.Width;
imheight = templateinfo.Height;

CellInd = find(cellpath{lasttp}(:,1)~=-1 & cellpath{lasttp}(:,2)~=-1)';
nucmask = cell(size(cellpath{lasttp},1),1);
cytomask = cell(size(cellpath{lasttp},1),1);
cellmask = cell(size(cellpath{lasttp},1),1);

for cellNo = CellInd
    load(celltrackOUTfilename,'cellpath','sisterList','bg');
    [pathLogic c_cellpath] = check_path(cellpath,sisterList,cellNo,firsttp,lasttp,imheight,imwidth);
    
    if pathLogic
        load(maskOUTfilename,'selected_cells');
        clc;display(['r' num2str(row) 'c' num2str(col) 'f' num2str(field) ',Processing cell:' num2str(cellNo)]);
        selected_cells = [selected_cells;cellNo];
        template = zeros(2*cellsize+1,2*cellsize+1,lasttp);
        
        newnucmask = zeros(size(template,1),size(template,2),lasttp);
        newcellmask = zeros(size(template,1),size(template,2),lasttp);
        newcytomask = zeros(size(template,1),size(template,2),lasttp);
        
        for tp = firsttp:lasttp
            if c_cellpath(tp,1)>0 && c_cellpath(tp,2)>0
                xL=max(c_cellpath(tp,1)-cellsize,1);
                xR=min(c_cellpath(tp,1)+cellsize,imwidth);
                yL=max(c_cellpath(tp,2)-cellsize,1);
                yR=min(c_cellpath(tp,2)+cellsize,imheight);
                imlocation = [row col field plane templateCH -1];
                
                if xR-xL == cellsize*2
                    borderX = 1:(cellsize*2+1);
                elseif xR == imwidth
                    shiftX = imwidth-xL;
                    borderX = 1:(shiftX+1);
                elseif xL == 1
                    shiftX = cellsize*2+1-xR;
                    borderX = (xL+shiftX):(cellsize*2+1);
                end
                
                if yR-yL == cellsize*2
                    borderY = 1:(cellsize*2+1);
                elseif yR == imheight
                    shiftY = imheight-yL;
                    borderY = 1:(shiftY+1);
                elseif yL == 1
                    shiftY = cellsize*2+1-yR;
                    borderY = (yL+shiftY):(cellsize*2+1);
                end
                
                template(borderY,borderX,tp) = loadimage(filetype,fileformat,imlocation,tp,channelnames,{[yL yR], [xL xR]},[]);
                
                [x1,y1,BW] = templateToCentroid(template(:,:,tp));
                
                if ~isempty(find(BW==1,1))
                    if c_cellpath(tp,1) > 0 && c_cellpath(tp,2) > 0

                        
                        newnucmask(:,:,tp) = BW;
                        newcellmask(:,:,tp) = bwmorph(BW,'dilate',celldilatesize);
                        newcytomask(:,:,tp) = bwmorph(newcellmask(:,:,tp),'erode',1) & ~newnucmask(:,:,tp);
                    end
                end
            end
        end
        nucName = ['nucmask_cell' num2str(cellNo)];
        cytoName = ['cytomask_cell' num2str(cellNo)];
        cellName = ['cellmask_cell' num2str(cellNo)];
        
        eval([nucName '=newnucmask;']);
        eval([cytoName '=newcytomask;']);
        eval([cellName '=newcellmask;']);
        
        save(maskOUTfilename,'selected_cells','-append');
        save(maskOUTfilename,nucName,'-append');
        save(maskOUTfilename,cytoName,'-append');
        save(maskOUTfilename,cellName,'-append');
        
        clear(nucName);
        clear(cytoName);
        clear(cellName);
    end
end

cd(currentF);

function [x y BW] = templateToCentroid(M)
BWc = zeros(size(M));
for i=1.2:0.6:3
    edgedIm = edge(M,'canny',0,i);
    BW = imfill(edgedIm,'holes');
    BW = bwmorph(BW,'open',1);
    BW = bwselect(BW,round(size(M,2)/2),round(size(M,1)/2));
    BWc = BWc | BW;
end
%figure(1),imshow(BWc)
S  = regionprops(BWc, 'centroid');
BW = zeros(size(M));
if isempty(find(BWc==0, 1)) || isempty(find(BWc==1, 1))
    x = round(size(M,2)/2);
    y = round(size(M,1)/2);
else
    x = round(S.Centroid(1));
    y = round(S.Centroid(2));
    [Ys,Xs] = find(BWc==1);
    xshift = x-round(size(M,2)/2);
    yshift = y-round(size(M,1)/2);
    if isempty(find( (Ys-yshift)<=0 | (Ys-yshift)>size(M,1) | (Xs-xshift)<=0 | (Xs-xshift)>size(M,2),1))
        for i=1:length(Ys)
            BW(Ys(i)-yshift,Xs(i)-xshift) = 1;
        end
    else
        BW = BWc;
    end
end


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

function [pathLogic new_cellpath] = check_path(cellpath,sisterList,cellNo,firsttp,lasttp,imheight,imwidth)
% check if this path is good
%#1 if acquiring and being inferior sister, combine trace of self with
%prior-sister
% has sister?
new_cellpath = zeros(lasttp,2);

if  sisterList{end}(cellNo,1) ~= -1
    sis1No = sisterList{lasttp}(cellNo,1);
    
    if cellpath{firsttp}(cellNo,1) == -1
        mainInd = sis1No;
        subInd =  cellNo;
        subData = zeros(1,lasttp);
        for t=firsttp:lasttp
            subData(t) = sisterList{t}(subInd,1);
        end
        t_break = find(subData(firsttp:lasttp)>0,1,'first');
        
        for t=firsttp:lasttp
            if t<t_break
                new_cellpath(t,:) = cellpath{t}(mainInd,:);
            else
                new_cellpath(t,:) = cellpath{t}(subInd,:);
            end
        end
    else
        for t=firsttp:lasttp
            new_cellpath(t,:) = cellpath{t}(cellNo,:);
        end
    end
else
    for t=firsttp:lasttp
        new_cellpath(t,:) = cellpath{t}(cellNo,:);
    end
end

deathInd = find(new_cellpath(:,1)==-2,1,'first');

if isempty(deathInd)
    temp_cellpath = new_cellpath;
else
    temp_cellpath = new_cellpath(1:(deathInd-1),:);
end

% Containing large x-y drift in comparison to background points ?

% Is there any point that lies outside the image area.
if isempty(  find(temp_cellpath(:,1)>imwidth | temp_cellpath(:,1)<0 | temp_cellpath(:,2)>imheight | temp_cellpath(:,2)<0,1)  )

    diff_X = diff(temp_cellpath(:,1));
    diff_Y = diff(temp_cellpath(:,2));
    med_XY = median([diff_X;diff_Y]);
    std_XY = std([diff_X;diff_Y]);
    
    if isempty(  find(diff_X > (med_XY+30*std_XY) | diff_Y > (med_XY+30*std_XY) ,1)  )
        pathLogic = 1;
    else
        pathLogic = 0;
    end
    
else
    display('Index out of bound');
    pathLogic = 0;
end



function outputim = loadimage(filetype,fileformat,imlocation,tp,channelnames,pixelsreg,displayGate)
row = imlocation(1);
col = imlocation(2);
field = imlocation(3);
plane = imlocation(4);
channel1 = imlocation(5);
channel2 = imlocation(6);
totalCH = length(channelnames);

if channel2 == -1
    switch filetype
        case 1
            filename = sprintf(fileformat,row,col,field,plane,channel1,tp);
            outputim = imread(filename,'PixelRegion',pixelsreg); 
        case 2
            outputim = imread(fileformat,'Index',totalCH*(tp-1)+channel1,'PixelRegion',pixelsreg);
        case 3
            filename = sprintf(fileformat,channelnames{channel1},tp);
            outputim = imread(filename,'PixelRegion',pixelsreg);
    end
    if isempty(displayGate)
        outputim = mat2gray(im2double(outputim));
    else
        ThresL = displayGate(1);
        ThresH = displayGate(2);
        outputim = mat2gray(im2double(outputim),[ThresL ThresH]);
    end
else
    switch filetype
        case 1
            filename = sprintf(fileformat,row,col,field,plane,channel1,tp);
            nominim = imread(filename,'PixelRegion',pixelsreg);
            filename = sprintf(fileformat,row,col,field,plane,channel2,tp);
            denomin = imread(filename,'PixelRegion',pixelsreg);
        case 2
            nominim = imread(fileformat,'Index',totalCH*(tp-1)+channel1,'PixelRegion',pixelsreg);
            denomin = imread(fileformat,'Index',totalCH*(tp-1)+channel2,'PixelRegion',pixelsreg);
        case 3
            filename = sprintf(fileformat,channelnames{channel1},tp);
            nominim = imread(filename,'PixelRegion',pixelsreg);
            filename = sprintf(fileformat,channelnames{channel2},tp);
            denomin = imread(filename,'PixelRegion',pixelsreg);
    end
    if isempty(displayGate)
        outputim = mat2gray(im2double(nominim+2^9)./im2double(denomin+2^9));
    else
        ThresL = displayGate(1);
        ThresH = displayGate(2);
        outputim = mat2gray(im2double(nominim+2^9)./im2double(denomin+2^9),[ThresL ThresH]);
    end
end

function iminfo = loadimageinfo(filetype,fileformat,imlocation,tp,channelnames,pixelsreg)
row = imlocation(1);
col = imlocation(2);
field = imlocation(3);
plane = imlocation(4);
channel = imlocation(5);
totalCH = length(channelnames);

switch filetype
    case 1
        filename = sprintf(fileformat,row,col,field,plane,channel,tp);
        iminfo = imfinfo(filename);
    case 2
        iminfo = imfinfo(fileformat,'Index',totalCH*(tp-1)+channel);
    case 3
        filename = sprintf(fileformat,channelnames{channel},tp);
        iminfo = imfinfo(filename);
end


