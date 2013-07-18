

% Define parameters related to the process---------
clear all;
signalshift = 0.01;
bgsubstractlogic = 0; % 
illumcorlogic = 0; % Algorithmic illumination correction by high-pass filter
framshift_logic = 0;
ImageIndex = 2; % 1=nomin/denomin, 2=templateCH, 3=nomin,4=denomin
intensityrange = [0.003357 0.23]; %#1 [0.0032807 0.008957] #2 [0 0.0084382] #3[0.0083314 0.22164]
displaygate = [0.85 1.6]; % For FRET Only
filterParam = [2 2]; 
cellsize = 15;
timestamplogic = 1; % 1 = frame no, 2 = actual time
celllocationlogic = 0; % 1 = show location of tracked cells, 0 = only image
save videoparameters;
clear all;
%-------------------------------------------------
% Define information about input images-----------
templateCH = 3; %1=foxO, 2=mCherry, 3=dye
nominCH = 2;
denominCH = 3;
sourcefolder = '/files/ImStor/sorger/data/Operetta/Bernhard/130520_184A1_CellTrackerDye_Stimulated/130520_CellTracker_Stimulated[261]/130520_CellTracker_Violet_Stimulation[516]/2013-05-20T210004Z[516]';
%------------------------------------------------


jobmgr = findResource('scheduler', 'type', 'lsf');
jobmgr.ClusterMatlabRoot = '/opt/matlab';
jobmgr.SubmitArguments = '-q short -W 12:00 -R "rusage[matlab_dc_lic=1]"';
job = jobmgr.createJob();

tps = [1 72];
sites = 1;

for row = 2:7
    for col = 1:12
        for field = 1
            
            fileformat = ['%03.0f%03.0f-%u-%03.0f001%03.0f.tif'];
            plane = 1;
            job.createTask(@GenMov_commandline, 0, ...
                {1,sourcefolder, row, col,field,plane,templateCH,nominCH,denominCH, tps,fileformat,[]});
            
        end
    end
end

job.submit();
