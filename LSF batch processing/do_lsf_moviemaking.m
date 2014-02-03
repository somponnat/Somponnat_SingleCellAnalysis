

% Define parameters related to the process---------
clear all;
signalshift = 0.0;
bgsubstractlogic = 0; % 
illumcorlogic = 0; % Algorithmic illumination correction by high-pass filter
framshift_logic = 0;
ImageIndex = 2; % 1=nomin/denomin, 2=templateCH, 3=nomin,4=denomin
intensityrange = [0.0053559 0.016]; %#1 [0.0032807 0.008957] #2 [0 0.0084382] #3[0.0083314 0.22164]
displaygate = [0.85 1.6]; % For FRET Only
filterParam = [2 2]; 
cellsize = 15;
timestamplogic = 2; % 1 = frame no, 2 = actual time
celllocationlogic = 0; % 1 = show location of tracked cells, 0 = only image
timestep_min = 5;
timestep_sec = 0;
save videoparameters;
clear all;
%-------------------------------------------------
% Define information about input images-----------
templateCH = 2; %1=foxO, 2=mCherry, 3=brightfield
nominCH = 1;
denominCH = 2;
sourcefolder = 'Q:\sorger\data\Operetta\Pat\12-10-2013 MCF10A FOXO3a construct timelapse with varying MEKi concentration\12-10-2013 MCF10A with varying MEKi concentration[730]\Post-drug-treatment[2291]\2013-12-10 175607 -0500[2301]';
%------------------------------------------------

tps = [1 200];
movieInd = [];
ind = 1;
for row = 2:4
    for col = 2:11
        for field = 1:3
            
            movieInd(ind,:) = [row col field];
            ind = ind+1;
        end
    end
end

if matlabpool('size') == 0
  matlabpool open;
end


parfor i = 1:length(movieInd)
    row = movieInd(i,1);
    col = movieInd(i,2);
    field = movieInd(i,3);
    fileformat = ['%03.0f%03.0f-%u-%03.0f001%03.0f.tif'];
    plane = 1;
    GenMov_commandline(1,sourcefolder,row,col,field,plane,templateCH,nominCH,denominCH,tps,fileformat,[]);
    
end


