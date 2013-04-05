jobmgr = findResource('scheduler', 'type', 'lsf');
jobmgr.ClusterMatlabRoot = '/opt/matlab';
jobmgr.SubmitArguments = '-q short -W 6:00 -R "rusage[matlab_dc_lic=1]"';

currentfolder = pwd;
targetfolder = '/files/ImStor/sorger/data/Operetta/Pat/12-10-2012/184A1/EGF20ngPerML/184a1-treatment__2012-12-10T14_04_00-Measurement1/Images';
rows = 5;
cols = 2:10;
fields = 1;
plane = 1;
channels = 2;
tps = [1 46];
job = jobmgr.createJob();
for row = rows
    for col = cols
        for field = fields
            for channel = channels
                cd(targetfolder);
                celltrackfile = ['celltrackOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(channel)];
                if exist([celltrackfile '.mat'],'file') == 2
                    
                    cd(currentfolder);
                    job.createTask(@CellTracking_commandline, 0, ...
                        {targetfolder, row, col, field, plane, channel, tps});
                else
                    cd(currentfolder);
                end
                
                
            end
        end
    end
end
job.submit();
job.wait();
job.destroy();
