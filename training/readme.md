Welcome to SUGAR training folder. Here are some useful scripts to prepare training dataset and the main scripts for training. 

We used the fk program to calculate the synthetic waveforms. For trial runs, use 10 earthquakes instead of 1000 to save time and debug. 

**Step 1 : run "to_fk_random_big.py"**

This will generate many outputs in the folder "random_event_allst". The initial generation of event locations and magnitudes have been commented. You need to modify the code. Those csv files numbered 0 to 999 are 1000 input earthquakes information with different location and magnitude for fk calculation. 

**Step 2 : Create a velocity model for fk**

The velocity model used by fk is a layered model. In each line, the values are 
thickness vs vp [rho Qs Qp]
Details are provided by the fk program. A used Chinese webpage is the Seisman note "https://blog.seisman.info/fk-notes/". 
The velocity model used in this study is "nz22vel", which was averaged from New Zealand wide 3D velocity model version 2.2. 

**Step 3 : run "GenSynEQ1000.py"**

This step will generate synthetic waveforms for the 1000 earthquakes. The output will be in the folder "fk_outputs". Create the folder before running. This script takes advantage of paralell computing. You need to set proper number of processes according to your computer capacity. 

**Step 4 : run "read_synthetic1000.py"**

The direct output of fk is in Radius and Transform directions. This step will change those to ZNE directions and save the files to the folder "fk_zne". Create the folder before running. 

**Step 5 : run "train500allst.py" and other "train##allst.py" files**

This step creates continous waveforms based on individual earthquake waveforms. The number 500 means there are about 500 earhtquakes per hour. Each file generate several hours of data in a given date, e.g., 1970-02-01, 1970-02-02, etc. We suggest to store the continuous waveforms in the folders formatted as "somewhere/data1970/02/01clean/". Each day has a different earthquake rate. Remember these are clean waveforms and we will add noise later. The input and output directory should be changed according to your folder structure. 

**Step 6 : run "train500noise.py" and others**

This step will add noise to the continuous waveform. In the training stage, we want to add noise in a range. So set the lower and upper boundaries and the program will give different hours in that day different noise levels. The waveforms are ground velocity in cm/s. So, in order to figure out a realistic noise level or level range, you need to see some real seismograms in a few stations, remove instrument response and filter to the ideal frequency band, and do the measurement. This does not need to be very accurate, because (a) we are using a range, and (b) we need to try a few different noise level ranges anyways to obtain the optimal trained model. We suggest to give the lower boundary much lower than the real noise level and the upper boundary slightly higher than the real one. The AI model should learn more clear samples than ambiguous samples. A useful way to see whether the noise level makes sense is to do the source scanning on some known earthquakes with real and synthetic data. Earthquakes with similar location and similar magnitude should result in similar brightness peaks in the brightness video. 

**Step 7 : run "toAIbrmaps_mul.py"**

This step calculates continuous brightness videos with those continuous waveforms. Note that you need to set how many hours in each day of the continuous waveforms. In the original application, some days have 20 hours and some days have 10 hours. You could set up everything around line 220 of this script. In this way, this one script can run through multiple days of data. In each hour, this script use paralell computing. Note that the continuous brightness video generated from this step is not yet the input training dataset. 

**Step 8 : Prepare synthetic earthquake catalogs to label the training dataset**

This can be done by running "eq1000catalog.py" for each day. The current example is for 1970-01-04, which contains 400 earthquakes. This was set in step 5. You need to change the total number of earthquakes according to what you set in step 5 for each day. After running this for all dates, each day should have a catalog stored in the clean data folder. 

**Step 9 : run "ssa2mat_smooth60_new"**

This step will cut out 60 by 60 by 60 blocks from the continuous brightness video and save corresponding 20 by 20 by 20 label blocks according to the input catalogs. We suggest to keep the folder structure as shown in this example. 

**Step 10 : run "delete0pad20to60.m"**

This step will delete empty samples and pad 0 to label blocks in order to get the same size of the training block. 

**Step 11 : Randomly put 10% data to validation set**

This can be achieved by running the following lines:

First enter the training data sample folder "imagesTr". 
```
ls | shuf -n 19000 | xargs -i mv {} ../imagesVal/
cd ../imagesVal/
ls *>filelist
mv filelist ../labelsTr/
cd ../labelsTr/
xargs -a filelist mv -t ../labelsVal/
rm filelist 
```
Then you could use the command "ls | wc -l" to see how many files are there in each folder. 

**Step 12 : run "eqtrain_linear_scaling.m"**

See settings inside the script. In this training process, you may see loss going done all the time, but overfitting may have already happened. In this study, we trained the U-Net 20 epochs to obtain the best result.

**Step 13 : Evaluate the trained model**

A testing dataset is required for model evaluation. You may used previous steps and scripts to generate a different set of continuous waveform, with some noise. We recommand having different earthquake rates in different dates, like training data, and in each day just one hour is okay. In testing mode, we don't need to cut the continuous brightness videos into blocks. Just save it as ".mat" format. The lines are commented at the moment in "toAIbrmaps_mul.py". 

The saved .mat brightness videos in one folder will be the inputs of "predictsmooth.m". This script will use a trained model to transfer the continuous brightness videos to continuous score videos. Note that the code needs to be changed with correct model name after training. Only the final saved model can be directly used in this script. If you want to use any epoch saved before the final model, you need to run "continuetrainsmooth.m" to format the checkpoint into a final model. 

After having those score videos, you could run "ai_catalog_f1.m" for each testing hour. This script will set a range of score threshold and output corresponding precision, recall and f1 values to a csv file. Note that the locations here are not finalized because we have not gone through phase picking and location refinement steps. But these solutions can provide a good estimates of the model performance. 

**Step 14 : Evaluate with real data**

You should have a few hours of real data ready and test the trained model with those real data. In that case, the entire workflow can be run and finalized locations and magnitudes can be obtained. The performance can be evaluated using the real data test. If the output catalog seems to contain too many detections and include many questionable detections, this means that the training dataset may have had too high noise level. In contrast, if the output catalog seems to miss too many events, this means that the training dataset may have had too low noise level. According to the performance in the real data test, you may go back to adjust the noise level range of the training dataset. Repeat the training a few times and select an optimal model. 






