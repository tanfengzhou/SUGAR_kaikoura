from obspy import read
import os
from obspy import Stream
from obspy import UTCDateTime
from pathlib import Path
import numpy as np
from EQTransformer.utils.hdf5_maker import preprocessor
from EQTransformer.core.predictor import predictor
import pickle
import pandas as pd
import shutil

def filter_times(pstation, pprob, sstation, sprob, within):
    all_times = set(pstation + sstation)
    for ptime, pprobability in zip(pstation, pprob):
        index = pstation.index(ptime)
        popped = 0
        for stime, sprobability in zip(sstation, sprob):
            if abs(stime - ptime) <= within and sprobability > pprobability:
                pstation.pop(index)
                pprob.pop(index)
                popped = 1
                break  
        
        if popped == 1:
            continue

        for other_time, other_probability in zip(pstation, pprob):
            if abs(other_time - ptime) <= within and other_time != ptime and other_probability > pprobability:
                pstation.pop(index)
                pprob.pop(index)
                break

    for stime, sprobability in zip(sstation, sprob):
        index = sstation.index(stime)
        popped = 0
        for ptime, pprobability in zip(pstation, pprob):
            if abs(ptime - stime) <= within and pprobability > sprobability:
                sstation.pop(index)
                sprob.pop(index)
                popped = 1
                break

        if popped == 1:
            continue

        for other_time, other_probability in zip(sstation, sprob):
            if abs(other_time - stime) <= within and other_time != stime and other_probability > sprobability:
                sstation.pop(index)
                sprob.pop(index)
                break

    return pstation, pprob, sstation, sprob


def main():

    time_begin=dict([
        ('year','2016'),
        ('month','11'),    # current version supports within 1 month only
        ('day','13'),
        ('hour','13'),
        ('min','00'),
        ('second', '00')
    ])

    time_end=dict([
        ('year','2016'),
        ('month','11'),
        ('day','13'),
        ('hour','14'),
        ('min','00'),     # keep it to the whole hour for the moment
        ('second', '00')
    ])

    process_unit = 60       # min
    overlap = 5             # min

    datadir = '/mnt/readynas6/ftan/2020location_method_leo2p/data2016all/'
    outputdir = './test1313eqt/'
    tempmseeddir = './temp/'

    totalday = int(float(time_end['day']) - float(time_begin['day']) + 1)
    for date in range(0, totalday):

        day = int(float(time_begin['day']) + date)
        if day < 10:
            day='0' + str(day)
        else:
            day=str(day)

        filenames = os.listdir(datadir + time_begin['month'] + '/' + day + '/')
        #print(filenames)

        if day == time_end['day']:
            duration = int(float(time_end['hour']))
        elif day == time_begin['day']:
            duration = 24 - int(float(time_begin['hour']))
        else:
            duration = 24

        if day == time_end['day'] and day == time_begin['day']:
            duration = int(float(time_end['hour'])) - int(float(time_begin['hour']))

        for runhour in range(0,duration):

            
            os.mkdir(tempmseeddir)

            if day == time_begin['day']:
                runh = runhour
            else:
                runh = 24 - int(float(time_begin['hour'])) + (date-1)*24 + runhour

            t1 = UTCDateTime(time_begin['year'] + '-' + time_begin['month'] + '-' + str(
                int(float(time_begin['day']) + np.floor(float(time_begin['hour']) + runh) / 24)) + 'T' + str(
                int(float(time_begin['hour']) + runh) % 24) + ':00:00')
            
            
            with open(outputdir + day + 'stations.p', 'rb') as f:
                comein = pickle.load(f)
            [stlas, stlos, stzs, stnms] = comein

            #print(stnms)

            for i in filenames:
                try:
                    channels = read(datadir + time_begin['month'] + '/' + day + '/' + i)
                    print('read in ' + i)
                except TypeError:
                    print('cannot read ' + i)
                    continue
                channels.merge()
                st = channels.split()
                st.detrend()
                st.merge(fill_value=0)

                if int(float(day) + 1) < 10:
                    dayp = '0' + str(int(float(day) + 1))
                else:
                    dayp = str(int(float(day) + 1))

                name1 = datadir + time_begin['month'] + '/' + dayp + '/' + time_begin['year'] + time_begin[
                    'month'] + dayp + '.' + st[0].stats.station + '.mseed'
                if Path(name1).is_file() == True:
                    try:
                        channelsnew = read(name1)
                        channelsnew.merge()
                        stnew = channelsnew.split()
                        stnew.detrend()
                        stnew.merge(fill_value=0)
                        combine = st + stnew
                        combine.merge(fill_value=0)
                        st = combine
                    except TypeError:
                        pass

                for j in st:
                    channel = j.stats.channel
                    stnm = j.stats.station
                    if stnm not in stnms:
                        print('station ' + stnm + ' not in the stations.p, skip it')
                        continue
                    else:
                        print('writing station ' + stnm)

                    try:
                        os.mkdir(tempmseeddir + stnm)
                    except FileExistsError:
                        pass

                    if j.stats.channel == 'BNZ' or j.stats.channel == 'HHZ' or j.stats.channel == 'EHZ':
                        newst = Stream([j])
                        newst.trim(t1, t1 + (process_unit + overlap) * 60,nearest_sample=False, pad=True, fill_value=0)
                        newst.write(tempmseeddir + stnm + '/NZ.' + stnm + '..' + channel + '__'
                                    + str(t1) + '__' + str(t1 + (process_unit + overlap) * 60) + '.mseed')

                    if (j.stats.channel == 'BN1' or j.stats.channel == 'BNE') or (
                            j.stats.channel == 'HH1' or j.stats.channel == 'HHE') or (
                            j.stats.channel == 'EH1' or j.stats.channel == 'EHE'):
                        newst = Stream([j])
                        newst.trim(t1, t1 + (process_unit + overlap) * 60,nearest_sample=False, pad=True, fill_value=0)
                        newst.write(tempmseeddir + stnm + '/NZ.' + stnm + '..' + channel + '__'
                                    + str(t1) + '__' + str(t1 + (process_unit + overlap) * 60) + '.mseed')

                    if (j.stats.channel == 'BN2' or j.stats.channel == 'BNN') or (
                            j.stats.channel == 'HH2' or j.stats.channel == 'HHN') or (
                            j.stats.channel == 'EH2' or j.stats.channel == 'EHN'):
                        newst = Stream([j])
                        newst.trim(t1, t1 + (process_unit + overlap) * 60,nearest_sample=False, pad=True, fill_value=0)
                        newst.write(tempmseeddir + stnm + '/NZ.' + stnm + '..' + channel + '__'
                                    + str(t1) + '__' + str(t1 + (process_unit + overlap) * 60) + '.mseed')

            json_basepath = os.path.join(os.getcwd(), "station_list_nz_real.json")

            ###########################################################
            # data preprocessing
            ###########################################################

            preprocessor(preproc_dir="temp2",
                         mseed_dir='temp',
                         stations_json=json_basepath,
                         overlap=0.7,
                         n_processor=100)

            ###########################################################
            # prediction
            ###########################################################

            predictor(input_dir='temp_processed_hdfs',
                      input_model='EqT_original_model.h5',
                      output_dir='temp_detections',
                      detection_threshold=0.3,
                      P_threshold=0.1,
                      S_threshold=0.1,
                      number_of_plots=1,
                      plot_mode='time',
                      #number_of_cpus=100,
                      output_probabilities=False)

            
            ###########################################################
            # save to sugar format
            ###########################################################

            finaldir = outputdir
            eqt_dir = './temp_detections/'
            sr = 50
            points = (process_unit + overlap) * 60 * sr

            with open(finaldir + day + str(int(float(time_begin['hour'])+runh)%24) + 'wood.p', 'rb') as f:
                comein = pickle.load(f)
            [vwood, h1wood, h2wood] = comein

            vwood = np.array(vwood)
            h1wood = np.array(h1wood)
            h2wood = np.array(h2wood)

            with open(finaldir + day + 'stations.p', 'rb') as f:
                comein = pickle.load(f)
            [stlas, stlos, stz, stnm] = comein

            nsta = len(stnm)

            onlyphaseP = np.zeros([nsta, points])
            onlyphaseh1 = np.zeros([nsta, points])
            onlyphaseh2 = np.zeros([nsta, points])
            Pmark = []
            Smark = []

            for i in stnm:
                stationdir = eqt_dir + i + '_outputs/'

                comein = pd.read_csv(stationdir + 'X_prediction_results.csv')
                pstation = []
                sstation = []
                pprob = []
                sprob = []

                for j in range(0, len(comein)):
                    row = comein.iloc[j]
                    ptime = row['p_arrival_time']
                    stime = row['s_arrival_time']
                    ppp = row['p_probability']
                    sss = row['s_probability']
                    try:
                        ptimehour = UTCDateTime(ptime) - t1
                        pstation.append(int(round(ptimehour * sr)))
                        pprob.append(float(ppp))
                    except TypeError:
                        pass

                    try:
                        stimehour = UTCDateTime(stime) - t1
                        sstation.append(int(round(stimehour * sr)))
                        sprob.append(float(sss))
                    except TypeError:
                        pass

                [pstation, pprob, sstation, sprob] = filter_times(pstation, pprob, sstation, sprob, within=sr * 1)

                Pmark.append(pstation)
                Smark.append(sstation)

            pam = []
            h1am = []
            h2am = []

            for i in range(nsta):
                station_am = []
                for j in range(len(Pmark[i])):
                    check_later = np.nonzero(onlyphaseP[i][Pmark[i][j] + sr: Pmark[i][j] + 5 * sr])
                    if len(check_later[0]) > 0:
                        nextphase = check_later[0][0] + sr
                        a = vwood[i][Pmark[i][j]: Pmark[i][j] + nextphase]
                    else:
                        a = vwood[i][Pmark[i][j]: Pmark[i][j] + 5 * sr]

                    am = max(abs(a))
                    onlyphaseP[i][Pmark[i][j]] = am
                    station_am.append(am)
                pam.append(station_am)

            for i in range(nsta):
                station_am = []
                for j in range(len(Smark[i])):
                    check_later = np.nonzero(onlyphaseh1[i][Smark[i][j] + sr: Smark[i][j] + 5 * sr])
                    if len(check_later[0]) > 0:
                        nextphase = check_later[0][0] + sr
                        a = h1wood[i][Smark[i][j]: Smark[i][j] + nextphase]
                    else:
                        a = h1wood[i][Smark[i][j]: Smark[i][j] + 5 * sr]

                    am = max(abs(a))
                    onlyphaseh1[i][Smark[i][j]] = am
                    station_am.append(am)
                h1am.append(station_am)

            for i in range(nsta):
                station_am = []
                for j in range(len(Smark[i])):
                    check_later = np.nonzero(onlyphaseh2[i][Smark[i][j] + sr: Smark[i][j] + 5 * sr])
                    if len(check_later[0]) > 0:
                        nextphase = check_later[0][0] + sr
                        a = h2wood[i][Smark[i][j]: Smark[i][j] + nextphase]
                    else:
                        a = h2wood[i][Smark[i][j]: Smark[i][j] + 5 * sr]

                    am = max(abs(a))
                    onlyphaseh2[i][Smark[i][j]] = am
                    station_am.append(am)
                h2am.append(station_am)

            onlyphaseS = np.maximum(onlyphaseh1, onlyphaseh2)
            sam = np.maximum(h1am, h2am)

            #pickle.dump([Pmark, Smark, onlyphaseP, onlyphaseS, pam, sam], open(outputdir + day + str(int(float(time_begin['hour'])+runh)%24) + 'eqtpicks.p', 'wb'))

            #shutil.rmtree('temp')
            #shutil.rmtree('temp2')
            #shutil.rmtree('temp_processed_hdfs')
            #shutil.rmtree('temp_detections')

if __name__ == "__main__":
    main()
