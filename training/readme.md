Welcome to SUGAR training folder. Here are some useful scripts to prepare training dataset and the main scripts for training. 

We used the fk program to calculate the synthetic waveforms. For trial runs, use 10 earthquakes instead of 1000 to save time and debug. 

Step 1 : run "to_fk_random_big.py". 

This will generate many outputs in the folder "random_event_allst". The initial generation of event locations and magnitudes have been commented. You need to modify the code. Those csv files numbered 0 to 999 are 1000 input earthquakes information with different location and magnitude for fk calculation. 

Step 2 : Create a velocity model for fk

The velocity model used by fk is a layered model. In each line, the values are 
thickness vs vp [rho Qs Qp]
Details are provided by the fk program. A used Chinese webpage is the Seisman note "https://blog.seisman.info/fk-notes/". 
The velocity model used in this study is "nz22vel", which was averaged from New Zealand wide 3D velocity model version 2.2. 

Step 3 : run "GenSynEQ1000.py"

This step will generate synthetic waveforms for the 1000 earthquakes. The output will be in the folder "fk_outputs". Create the folder before running. This script takes advantage of paralell computing. You need to set proper number of processes according to your computer capacity. 

Step 4 : run "read_synthetic1000.py"

The direct output of fk is in Radius and Transform directions. This step will change those to ZNE directions and save the files to the folder "fk_zne". Create the folder before running. 

Step 5 : run "train500allst.py" and other "train##allst.py" files.

This step creates continous waveforms based on individual earthquake waveforms. The number 500 means there are about 500 earhtquakes per hour. Each file generate several hours of data in a given date, e.g., 1970-02-01, 1970-02-02, etc. We suggest to store the continuous waveforms in the folders formatted as "somewhere/data1970/02/01clean/". Each day has a different earthquake rate. Remember these are clean waveforms and we will add noise later. The input and output directory should be changed according to your folder structure. 

Step 6 : run "train500noise.py" and others. 

This step will add noise to the continuous waveform. In the training stage, we want to add noise in a range. So set the lower and upper boundaries and the program will give different hours in that day different noise levels. The waveforms are ground velocity in cm/s. So, in order to figure out a realistic noise level or level range, you need to see some real seismograms in a few stations, remove instrument response and filter to the ideal frequency band, and do the measurement. This does not need to be very accurate, because (a) we are using a range, and (b) we need to try a few different noise level ranges anyways to obtain the optimal trained model. We suggest to give the lower boundary much lower than the real noise level and the upper boundary slightly higher than the real one. The AI model should learn more clear samples than ambiguous samples. A useful way to see whether the noise level makes sense is to do the source scanning on some known earthquakes with real and synthetic data. Earthquakes with similar location and similar magnitude should result in similar brightness peaks in the brightness video. 

Step 7 : run "toAIbrmaps_mul.py"

This step calculates continuous brightness videos with those continuous waveforms. Note that you need to set how many hours in each day of the continuous waveforms. In the original application, some days have 20 hours and some days have 10 hours. You could set up everything around line 220 of this script. In this way, this one script can run through multiple days of data. In each hour, this script use paralell computing. Note that the continuous brightness video generated from this step is not yet the input training dataset. 

Step 8 : Prepare synthetic earthquake catalogs to label the training dataset. 

This can be done by running "eq1000catalog.py" for each day. The current example is for 1970-01-04, which contains 400 earthquakes. This was set in step 5. You need to change the total number of earthquakes according to what you set in step 5 for each day. After running this for all dates, each day should have a catalog stored in the clean data folder. 

Step 9 : run "ssa2mat_smooth60_new"

This step will cut out 60*60*60 blocks from the continuous brightness video and save corresponding 20*20*20 label blocks according to the input catalogs. We suggest to keep the folder structure as shown in this example. 

Step 10 : run "delete0pad20to60.m" 

This step will delete empty samples and pad 0 to label blocks in order to get the same size of the training block. 

Step 11 : Randomly put 10% data to validation set. 

This can be achieved by running the following lines:

First enter the training data sample folder "imagesTr". 

  ls | shuf -n 19000 | xargs -i mv {} ../imagesVal/
  cd ../imagesVal/ 
  ls *>filelist
  mv filelist ../labelsTr/
  cd ../labelsTr/
  xargs -a filelist mv -t ../labelsVal/
  rm filelist 

Then you could use the command "ls | wc -l" to see how many files are there in each folder. 





