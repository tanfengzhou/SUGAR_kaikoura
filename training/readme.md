Welcome to SUGAR training folder. Here are some useful scripts to prepare training dataset and the main scripts for training. 

We used the fk program to calculate the synthetic waveforms. 

Step 1 : run "to_fk_random_big.py". 

This will generate many outputs in the folder "random_event_allst". The initial generation of event locations and magnitudes have been commented. You need to modify the code. Those csv files numbered 0 to 999 are 1000 input earthquakes information with different location and magnitude for fk calculation. 

Step 2 : Create a velocity model for fk

The velocity model used by fk is a layered model. In each line, the values are 
thickness vs vp [rho Qs Qp]
Details are provided by the fk program. A used Chinese webpage is the Seisman note "https://blog.seisman.info/fk-notes/". 
The velocity model used in this study is "nz22vel", which was averaged from New Zealand wide 3D velocity model version 2.2. 

Step 3 : run "GenSynEQ1000.py"

This step will generate synthetic waveforms for the 1000 earthquakes. The output will be in the folder "fk_outputs". This script takes advantage of paralell computing. You need to set proper number of processes according to your computer capacity. 
