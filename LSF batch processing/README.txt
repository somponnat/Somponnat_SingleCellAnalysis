1. Optionally, create a suitable <MATLAB> directory:

  % mkdir -p <MATLAB>

2. Copy files to the <MATLAB> directory as shown below:

  % cp -a /home/gfb2/share/ss240/CellTracking_commandline.m <MATLAB>
  % cp -a /home/gfb2/share/ss240/do_lsf.m <MATLAB>

3. Make <MATLAB> the current directory:

  % cd <MATLAB>
  
4. Start MATLAB as follows:

  % matlab -nodisplay -nodesktop -nosplash

5. Run do_lsf at the MATLAB command line:

  >> do_lsf

   If all goes well, after all the jobs finish (it should take 30-60 seconds),
   the screen should look like this (though with a different job number):

  >> do_lsf
  bkill 151288: Signal 127

6. Optionally, in a separate terminal window, run the following command:

  % watch bjobs

   ...to see the current status of the spawned jobs.  (You may want to execute
   this command *before* running do_lsf in the MATLAB session.)  At first, the
   output of `watch bjobs` will be only a couple of lines (refreshed every 2
   seconds):

  Every 2.0s: bjobs                                       Fri Nov 23 11:04:05 2012

  No unfinished job found


   ...but eventually the output will look like this:

  Every 2.0s: bjobs                                       Fri Nov 23 17:05:09 2012

  JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[21]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[6]   Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[24]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet002 Job14[12]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[8]   Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[19]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[7]   Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet002 Job14[14]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet002 Job14[13]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[5]   Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[4]   Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet002 Job14[2]   Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[17]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[1]   Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[18]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[9]   Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[20]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet002 Job14[3]   Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet002 Job14[10]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet002 Job14[11]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[15]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[16]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[22]  Nov 23 17:04
  151288  ss240   RUN   shared_2h  mezzanine.o clarinet001 Job14[23]  Nov 23 17:04

   To quit `watch bsubs`, just hit Ctrl-C.

