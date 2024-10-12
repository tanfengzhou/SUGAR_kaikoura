#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 10:01:11 2019
@author: Dawei Gao, Fengzhou Tan
Modified from: ArrayScriptDis_DGv2.py
Modified from: GenSyn.py
"""

# %%
import numpy as np
import math
import pandas as pd
import os, shutil
import time
import multiprocessing

def round_dec2(A):
    return float("{:.2f}".format(A))

def round_dec5(A):
    return float("{:.2f}".format(A))

Velmodel = 'nz22vel'
FilePath = './random_event_allst/'

EQrange = range(0, 1000)

def multi_fk(EQid):

    ResFolder = './fk_outputs/EQ{0}/'.format(EQid)
    if os.path.exists(ResFolder):
        os.system('rm -r -f {0}'.format(ResFolder))
    os.makedirs(ResFolder)

    stafolder = ResFolder + 'StaInfo/'
    os.makedirs(stafolder)

    StaFile = '{0}.csv'.format(EQid)
    df = pd.read_csv(FilePath+StaFile,names=['StaID_old', 'SRdis', 'Azimuth', 'Elevation','EQdep','Strike','Dip','Rake','SourceDuration','EQmag'])

    df.sort_values(by=['SRdis'],inplace=True)
    df.reset_index(drop=True,inplace=True)
    StaID       = df.index.tolist()
    df['StaID'] = StaID    # renew station ID

    vmodel = Velmodel
    npts   = '8192'
#    npts   = '4096'
    deltaT = '0.025'

    mag    = '{0}'.format(df['EQmag'] .tolist()[0])
    strike = '{0}'.format(df['Strike'].tolist()[0])
    dip    = '{0}'.format(df['Dip']   .tolist()[0])
    rake   = '{0}'.format(df['Rake']  .tolist()[0])

    duration = [round_dec5(df['SourceDuration'].tolist()[0])]
    risetime = 0.5

    EQdep = '{0}'.format(df['EQdep'].tolist()[0])     # EQ depth

    refdis = '1000.0'  # reference Distance

    DisAzifloat = 3 # precision for distance and azimuth for rest stations
    n_dis       = 80 # number of distances each line for Green's calculation

    Flag0 = 0 # waveform: velocity     [default]
    #Flag0 = 1 # waveform: displacement

    SRdis       = df['SRdis'].tolist()
    Azimuth     = df['Azimuth'].tolist()

    SRdis   = [str(round(x,DisAzifloat)) for x in SRdis  ]
    Azimuth = [str(round(x,DisAzifloat)) for x in Azimuth]

    df['SRdis']   = SRdis
    df['Azimuth'] = Azimuth
    df['SourceDuration'] = df['SourceDuration'].apply(round_dec5)
    df.to_csv(ResFolder+'StaInfo.csv',index=False)

    NumSta = len(StaID)

    #==============================================================================
    # Write file
    #==============================================================================

    if Flag0 == 1:
        text_file = open('./fk_outputs/EQ' + str(EQid) + '/runSynEQdis.csh', "w")
    else:
        text_file = open('./fk_outputs/EQ' + str(EQid) + '/runSynEQvel.csh', "w")

    # parameters
    text_file.write("#!/bin/csh\n")
    text_file.write("\n")
    text_file.write("#echo ''\n")
    text_file.write("#echo 'Generating waveforms... There are {0} stations in total.'\n".format(NumSta))
    text_file.write("\n")
    text_file.write("set vmodel = {0}\n".format(vmodel))
    text_file.write("set depth  = {0}\n".format(EQdep))
    text_file.write("set npts   = {0}\n".format(npts))
    text_file.write("set deltaT = {0}\n".format(deltaT))
    text_file.write("\n")

    text_file.write("mkdir ./WFraw\n")

    text_file.write("\n")
    text_file.write("#echo ''\n")
    text_file.write("#date\n")
    text_file.write("#echo ''\n")
    text_file.write("\n")
    text_file.write("\n")
    text_file.write("#---------------------------------- calculate green's funtions ----------------------------------#\n")

    #-- green's function
    n_job = math.ceil(len(SRdis)/n_dis)
    for j in range(n_job):
        if j < n_job-1:
            tempdis = SRdis[j*n_dis:(j+1)*n_dis]
            line = 'fk.pl -M${vmodel}/${depth} -N${npts}/${deltaT}'
            for k in range(len(tempdis)):
                line = line + ' '+ tempdis[k]
            line = line + ' {0}\n'.format(refdis)
            text_file.write(line)
            text_file.write("\n")
        else:
            tempdis = SRdis[j*n_dis:]
            line = 'fk.pl -M${vmodel}/${depth} -N${npts}/${deltaT}'
            for k in range(len(tempdis)):
                line = line + ' '+ tempdis[k]
            line = line + ' {0}\n'.format(refdis)
            text_file.write(line)
            text_file.write("\n")
    text_file.write('mv junk.* ./{0}/\n'.format('${vmodel}_${depth}'))
    text_file.write("\n")

    #-- waveform
    text_file.write("#----------------------------------     calculate waveforms    ----------------------------------#\n")

    for i in range(NumSta):
        staID          = StaID[i]
        distance       = SRdis[i]
        azimuth        = Azimuth[i]
        text_file.write("\n")
        text_file.write("#------ Station {0}:    distance = {1} km    azimuth = {2} deg\n".format(staID,distance,azimuth))
        text_file.write("\n")
        station_folder = "./WFraw/Sta{0}_Azi{1}".format(StaID[i],Azimuth[i])
        text_file.write("mkdir {0}\n".format(station_folder))
        text_file.write("\n")

        for j in range(len(duration)):
            text_file.write("#({0}) source duration: {1}\n".format(j,duration[j]))

            if Flag0 == 1:
                print('')
                print('Displacement: check syn parameters please!!!')
                print('')
                pass

            else:
                junk = '${vmodel}_${depth}'
                line = 'syn -M{0}/{1}/{2}/{3} -D{4}/{5} -A{6} -Osta{7}.z -G'.format(mag,strike,dip,rake,duration[j],risetime,azimuth,staID)+junk+'/{0}.grn.0\n'.format(distance)
                text_file.write(line)

            text_file.write("mkdir {0}/WF{1}\n".format(station_folder,duration[j]))
            text_file.write("mv sta* {0}/WF{1}/\n".format(station_folder,duration[j]))
            text_file.write("\n")

        #-- save station info.
        df1 = pd.DataFrame()
        df1['Sta{0}_Dis'.format(staID)] = [distance]
        df1['Sta{0}_Azi'.format(staID)] = [azimuth]
        df1.to_csv(stafolder + 'Sta{0}_Azi{1}.csv'.format(staID,azimuth),index=False)

    text_file.write("#echo ''\n")
    text_file.write("#echo 'Calculation is done!'\n")
    text_file.write("#date\n")
    text_file.close()

    if Flag0 == 1:
        os.system('chmod +x ./fk_outputs/EQ' + str(EQid) + '/runSynEQdis.csh')
    else:
        os.system('chmod +x ./fk_outputs/EQ' + str(EQid) + '/runSynEQvel.csh')

    # note: this version only deals with velocity seismogram
    #os.system('mv ./runSynEQvel.csh ./fk_outputs/EQ{0}/'.format(EQid))
    os.system('cp ./{0} ./fk_outputs/EQ{1}/'.format(Velmodel,EQid))

    # generate waveforms
    owd = os.getcwd()
#    print(owd)
    os.chdir('./fk_outputs/EQ{0}/'.format(EQid))
#    print(os.getcwd())
    os.system('./runSynEQvel.csh')
    os.chdir(owd)
#    print(os.getcwd())
    print('EQ' + str(EQid) + ' is done')
    time.sleep(2)


if __name__ == '__main__':

    pool = multiprocessing.Pool(100)

    start_time = time.time()

    pool.map(multi_fk, EQrange)

    print("Time elapsed: ", round_dec2((time.time() - start_time)/3600), "hours.")

