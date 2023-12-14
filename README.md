Source Untangler Guided by Artificial intelligence image Recognition (SUGAR) 

Step 1, 2, 2.5, 3, 4 should be run sequencially. 
Different versions are available for each step. For example, in step 1, either a 1D or a 3D velocity model can be used. The "3D" in the file name indicates that the script is using a 3D velocity model. Step 2.5 is phase picking via Earthquake Transformer (Mousavi et al., 2020, Nature Communications). If in step 3 you choose kurtosis picking technique, step 2.5 should be skipped. 

Example data are available in "data2016all". Full data should be downloaded directly through the GNS website and modified into the format shown here. 

3D travel time tables should be calculated before running the main code. Scripts are available in "calculate_3d_time". The tables are too big to be uploaded. A computer with >64G memory may be required for the 3D traveltime tables. 

Training scripts and data examples are available in "training". For training, brightness videos should be generated and then cut to $60\times60\times60$ blocks. Labels should correspond to the center $20\times20\times20$ areas, and then run "delete0pad20to60.m" to prepare the data for training. The main training code is "eqtrain_linear_scaling.m". Matlab and a GPU with >8G memory are required. The 3D UNet structure is also available in other packages in Python, which may be used. 

Please contact Fengzhou Tan via fengzhou.tan@gmail.com if you have any question. 
