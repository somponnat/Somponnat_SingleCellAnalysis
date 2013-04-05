clear all;

cellsize = 15;
signalshift = 2^9;
bgsubstractlogic = 0; % 
illumcorlogic = 0;
framshift_logic = 0;
ImageIndex = 2; % 1=nomin/denomin, 2=template, 3=nomin,4=denomin
filterParam = [2 2];

intensityrange = [0 2000];
displaygate = [0.97 1.5];
timestep_min =20; %minutes
timestep_sec = 0;%second
timestamplogic = 2; % 1 = frame no, 2 = actual time
celllocationlogic = 0; % 1 = show location of tracked cells, 0 = only image

save videoparameters;
clear all;


jobmgr = findResource('scheduler', 'type', 'lsf');
jobmgr.ClusterMatlabRoot = '/opt/matlab';
jobmgr.SubmitArguments = '-q short -W 12:00 -R "rusage[matlab_dc_lic=1]"';

currentfolder = pwd;

targetfolder = '/files/ImStor/sorger/data/Operetta/Bernhard/121222_Starvation_Experiment_Cell_Lines_and_Media_Conditions/122212_Starvation_Experiment_184A1_FoxO3a_Cherry_EKAREV_2__2012-12-22T16_04_14-Measurement1/Images';
rows = 2:7;
cols = [5 6 9 10];


fields = 1:2;
plane = 1;
channels = [3];
tps = [1 50];
fileformat = 'r%02.0fc%02.0ff%02.0fp%02.0frc%1.0f-ch1sk%ufk1fl1.tiff';


job = jobmgr.createJob();

for row = rows
    for col = cols
        for field = fields
            for channel = channels     
                celltrackfile = ['celltrackOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(channel)];
                videoname = ['myMov_r' num2str(row) 'c' num2str(col) 'f' num2str(field) 'ch' num2str(channel) '.avi'];
                
                job.createTask(@GenMov_commandline, 0, ...
                    {1,targetfolder, row, col,field,plane,channel,-1,-1, tps,fileformat,~});
            end
        end
    end
end

job.submit();


