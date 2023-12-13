import numpy as np
import pickle
import csv
import geodis
import obspy
from obspy import read
from pathlib import Path
import math
import copy
import ssa
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.io import savemat
from woodanderson import wa
from obspy.core.stream import Stream
from multiprocessing import Pool

def process_trace(args):
    inv, i, j, k, sr, points = args
    my_tr = Stream(traces=[i, j, k])
    print(k)
    try:
        my_tr.remove_response(inventory=inv, output="VEL")
        # my_tr.filter('bandpass', freqmin=0.5, freqmax=10, zerophase=True)
        my_tr.detrend(type='demean')

        delta = 1 / sr

        for i in range(0, len(my_tr)):
            my_tr[i] = wa(my_tr[i], delta)

        return (my_tr[0].data, my_tr[1].data, my_tr[2].data)

    except (KeyboardInterrupt, SystemExit):
        raise

    except:
        print('fail')
        return (np.ones(points), np.ones(points), np.ones(points))


def parallel_process_traces(xmlfiles, h1raw, h2raw, vraw, sr, points, num_processes):
    invs = xmlfiles
    pool = Pool(num_processes)
    results = pool.map(process_trace, zip(invs, h1raw, h2raw, vraw, [sr] * len(invs), [points] * len(invs)))
    pool.close()
    pool.join()
    h1wood, h2wood, vwood = zip(*results)
    return np.array(h1wood), np.array(h2wood), np.array(vwood)


def main():
    ######################################################################
    # parameters
    ######################################################################

    nametag=''  # tag for different experiments

    #############################################################
    # set study area with margin of 20
    #############################################################

    studyarea = [-44.2, 171.0, 3.72, 5.2]
    studydepth=[10,10]
    xygrid=4         #km
    depgrid=4        #km

    #############################################################
    # set study period
    #############################################################

    time_begin=dict([
        ('year','2016'),
        ('month','12'),    # current version supports within 1 month only
        ('day','01'),
        ('hour','00'),
        ('min','00'),
        ('second', '00')
    ])

    time_end=dict([
        ('year','2016'),
        ('month','12'),
        ('day','31'),
        ('hour','24'),
        ('min','00'),     # keep it to the whole hour for the moment
        ('second', '00')
    ])

    process_unit = 60       # min
    overlap = 5             # min

    #############################################################
    # set folder and file names
    #############################################################


    datadir = '/mnt/readynas6/ftan/2020location_method_leo2p/data2016all/'
    outputdir = '12_3d/'
    xmldir = '/mnt/readynas6/ftan/2020location_method_leo2p/station_xml_new/'
    station_info_file = 'station_master_new.csv'
    area = 'fewer'     # area code in station info file
    traveltime_tableP_name = './calculate_3d_time/ptraveltimes3d_0.2.p'
    traveltime_tableS_name = './calculate_3d_time/straveltimes3d_0.2.p'

    #############################################################
    # set scanning parameters
    #############################################################

    sr=50      # Hz
    scanlowf=5     # frequency band, Hz
    scanhighf=20
    root=3        # don't change
    win=6       # scanning window duration, seconds
    step=0.5      # seconds
    processes = 120

    ######################################################################
    # main program starts
    ######################################################################

    points=int(sr*(process_unit+overlap)*60)

    with open(traveltime_tableP_name, 'rb') as f:
        comein=pickle.load(f)
    ptraveltimes3d = comein

    with open(traveltime_tableS_name, 'rb') as f:
        comein=pickle.load(f)
    straveltimes3d = comein

    complete_stnms = []
    with open('stnm_fewer', 'r') as f:
        comein = f.readlines()
    for i in comein:
        complete_stnms.append(i[:-1])

    latgrid=xygrid/geodis.geodis([[studyarea[0]-0.5,studyarea[1]],[studyarea[0]+0.5,studyarea[1]]])     #degree
    longrid=xygrid/geodis.geodis([[studyarea[0],studyarea[1]-0.5],[studyarea[0],studyarea[1]+0.5]])     #degree

    stations=[]
    with open(station_info_file,'r', encoding='utf-8-sig') as f:
        comein = csv.reader(f, delimiter='|')
        line_count=0
        for row in comein:
            if line_count==0:
                title=row[0].split(',')
                line_count=line_count+1
            else:
                stations.append(row[0].split(','))
                line_count=line_count+1

    area_index=title.index('area')
    net_index=title.index('net')
    name_index=title.index('sta')
    z_index=title.index('z')
    e_index=title.index('e')
    n_index=title.index('n')
    lat_index=title.index('lat')
    lon_index=title.index('lon')
    elev_index=title.index('elev')

    station_use=[]
    for row in stations:
        if row[area_index]==area:
            station_use.append(row)

    lats=np.arange(studyarea[0],studyarea[0]+studyarea[2],latgrid)
    lons=np.arange(studyarea[1],studyarea[1]+studyarea[3],longrid)
    deps=np.arange(studydepth[0],studydepth[1]+depgrid,depgrid)
    studygrids=[]
    for i in lats:
        for j in lons:
            for k in deps:
                studygrids.append([i,j,k])

    print(len(studygrids))
    pickle.dump([studygrids],open(outputdir + 'studygrids.p','wb'))


    totalday = int(float(time_end['day']) - float(time_begin['day']) + 1)
    for date in range(0, totalday):

        day = int(float(time_begin['day']) + date)
        if day < 10:
            day='0' + str(day)
        else:
            day=str(day)

        stlas = []
        stlos = []
        stz = []
        stnm = []

        v=[]
        h1=[]
        h2=[]

        xmlfiles=[]

        for i in station_use:
            name1 = datadir + time_begin['month'] + '/' + day + '/' + time_begin['year'] + time_begin['month'] + day + '.' + i[name_index] + '.mseed'
            if Path(name1).is_file() == True:
                try:
                    channels=read(name1)
                except TypeError:
                    continue

                channels.merge()
                stnew = channels.split()
                stnew.detrend()
                stnew.merge(fill_value=0)

                zcount = 0
                ncount = 0
                ecount = 0
                t1 = UTCDateTime(time_begin['year'] + '-' + time_begin['month'] + '-' + day + 'T00:00:00')
                for j in stnew:
                    if (j.stats.channel == 'BNZ' or j.stats.channel == 'HHZ' or j.stats.channel == 'EHZ') and zcount == 0:
                        if j.stats.sampling_rate != sr:
                            j.resample(sr)
                        v.append(j)
                        print(j)
                        zcount = 1
                    if ((j.stats.channel == 'BN1' or j.stats.channel == 'BNE') or (j.stats.channel == 'HH1' or j.stats.channel == 'HHE') or (j.stats.channel == 'EH1' or j.stats.channel == 'EHE')) and ncount == 0:
                        if j.stats.sampling_rate != sr:
                            j.resample(sr)
                        h1.append(j)
                        print(j)
                        ncount = 1
                    if ((j.stats.channel == 'BN2' or j.stats.channel == 'BNN') or (j.stats.channel == 'HH2' or j.stats.channel == 'HHN') or (j.stats.channel == 'EH2' or j.stats.channel == 'EHN')) and ecount == 0:
                        if j.stats.sampling_rate != sr:
                            j.resample(sr)
                        h2.append(j)
                        print(j)
                        ecount = 1

                if zcount == 1 and ncount == 1 and ecount==1:
                    stlas.append(float(i[lat_index]))
                    stlos.append(float(i[lon_index]))
                    stz.append(float(i[elev_index]))
                    stnm.append(i[name_index])
                    inv = obspy.read_inventory(xmldir +i[name_index] +'.xml', format='STATIONXML')
                    xmlfiles.append(inv)

        ptraveltimes = np.zeros((len(studygrids), len(stnm)))
        straveltimes = np.zeros((len(studygrids), len(stnm)))
        for i in range(len(stnm)):
            for j in range(len(complete_stnms)):
                if stnm[i] == complete_stnms[j]:
                    ptraveltimes[:, i] = ptraveltimes3d[:, j]
                    straveltimes[:, i] = straveltimes3d[:, j]

        lat0=min(stlas)-0.5
        lat=max(stlas)+0.5
        lon0=min(stlos)-1
        lon=max(stlos)+1

        m = Basemap(llcrnrlon=lon0,llcrnrlat=lat0,urcrnrlon=lon,urcrnrlat=lat,projection='cass',lat_0=(lat0+lat)/2,lon_0=(lon0+lon)/2)

        m.drawcoastlines()
        m.drawparallels(np.arange(lat0, lat, (lat-lat0)/4), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(lon0, lon, (lon-lon0)/4), labels=[0, 0, 0, 1])
        m.drawmapboundary()

        x, y = m(stlos, stlas)
        m.scatter(x, y, 40, color='g', marker='^')
        #x, y = m(-129.48, 49.25)
        #m.scatter(x, y, 100, color='r', marker='*')

        plt.savefig(outputdir + day +'stations.pdf')

        pickle.dump([stlas, stlos, stz, stnm], open(outputdir + day + 'stations.p', 'wb'))

        #######################################################################

        nsta=len(v)

        if day == time_end['day']:
            duration = int(float(time_end['hour']))
        elif day == time_begin['day']:
            duration = 24 - int(float(time_begin['hour']))
        else:
            duration = 24

        if day == time_end['day'] and day == time_begin['day']:
            duration = int(float(time_end['hour'])) - int(float(time_begin['hour']))

        for runhour in range(0,duration):

            if day == time_begin['day']:
                runh = runhour
            else:
                runh = 24 - int(float(time_begin['hour'])) + (date-1)*24 + runhour

            if int(float(time_begin['hour'])+runh)%24 == 23:
                vborrow = []
                h1borrow = []
                h2borrow = []
                if int(float(day)+1) < 10:
                    dayp = '0' + str(int(float(day)+1))
                else:
                    dayp = str(int(float(day)+1))
                for i in station_use:
                    name1 = datadir + time_begin['month'] + '/' + dayp + '/' + time_begin['year'] + time_begin['month'] + dayp + '.' + i[name_index] + '.mseed'
                    if Path(name1).is_file() == True:
                        try:
                            channels=read(name1)
                        except TypeError:
                            continue
                        #channels.merge(method=1, fill_value='interpolate')
                        channels.merge()
                        stnew = channels.split()
                        stnew.detrend()
                        stnew.merge(fill_value=0)

                        zcount = 0
                        ncount = 0
                        ecount = 0
                        t1 = UTCDateTime(time_begin['year'] + '-' + time_begin['month'] + '-' + day + 'T00:00:00')
                        for j in stnew:
                            if (j.stats.channel == 'BNZ' or j.stats.channel == 'HHZ' or j.stats.channel == 'EHZ') and zcount == 0:
                                if j.stats.sampling_rate != sr:
                                    j.resample(sr)
                                vborrow.append(j)
                                print(j)
                                zcount = 1
                            if ((j.stats.channel == 'BN1' or j.stats.channel == 'BNE') or (j.stats.channel == 'HH1' or j.stats.channel == 'HHE') or (j.stats.channel == 'EH1' or j.stats.channel == 'EHE')) and ncount == 0:
                                if j.stats.sampling_rate != sr:
                                    j.resample(sr)
                                h1borrow.append(j)
                                print(j)
                                ncount = 1
                            if ((j.stats.channel == 'BN2' or j.stats.channel == 'BNN') or (j.stats.channel == 'HH2' or j.stats.channel == 'HHN') or (j.stats.channel == 'EH2' or j.stats.channel == 'EHN')) and ecount == 0:
                                if j.stats.sampling_rate != sr:
                                    j.resample(sr)
                                h2borrow.append(j)
                                print(j)
                                ecount = 1

            t1 = UTCDateTime(time_begin['year']+'-'+time_begin['month']+'-'+str(int(float(time_begin['day'])+np.floor(float(time_begin['hour'])+runh)/24))+'T'+str(int(float(time_begin['hour'])+runh)%24)+':00:00')

            vfilter = []
            h1filter = []
            h2filter = []

            vraw = []
            h1raw = []
            h2raw = []

            for i in v:
                ii = copy.deepcopy(i)
                if int(float(time_begin['hour'])+runh)%24 == 23:
                    for look in vborrow:
                        if look.stats.station == ii.stats.station:
                            st = Stream(traces=[ii, look])
                            st.merge(fill_value=0)
                            ii = st[0]
                            break
                try:
                    ii.trim(t1, t1 + (process_unit + overlap) * 60, nearest_sample=False, pad=True, fill_value=0)
                    ii.detrend()
                    ii.filter(type='bandpass', freqmin=scanlowf, freqmax=scanhighf, corners=3, zerophase=True)
                    vfilter.append(ii.data)
                except NotImplementedError:
                    vfilter.append(np.ones(points))

            for i in h1:
                ii = copy.deepcopy(i)
                if int(float(time_begin['hour']) + runh) % 24 == 23:
                    for look in h1borrow:
                        if look.stats.station == ii.stats.station:
                            st = Stream(traces=[ii, look])
                            st.merge(fill_value=0)
                            ii = st[0]
                            break
                try:
                    ii.trim(t1, t1 + (process_unit + overlap) * 60, nearest_sample=False, pad=True, fill_value=0)
                    ii.detrend()
                    ii.filter(type='bandpass', freqmin=scanlowf, freqmax=scanhighf, corners=3, zerophase=True)
                    h1filter.append(ii.data)
                except NotImplementedError:
                    h1filter.append(np.ones(points))

            for i in h2:
                ii = copy.deepcopy(i)
                if int(float(time_begin['hour']) + runh) % 24 == 23:
                    for look in h2borrow:
                        if look.stats.station == ii.stats.station:
                            st = Stream(traces=[ii, look])
                            st.merge(fill_value=0)
                            ii = st[0]
                            break
                try:
                    ii.trim(t1, t1 + (process_unit + overlap) * 60 ,nearest_sample=False, pad=True, fill_value=0)
                    ii.detrend()
                    ii.filter(type='bandpass', freqmin=scanlowf, freqmax=scanhighf, corners=3, zerophase=True)
                    h2filter.append(ii.data)
                except NotImplementedError:
                    h2filter.append(np.ones(points))


            for i in v:
                ii = copy.deepcopy(i)
                if int(float(time_begin['hour'])+runh)%24 == 23:
                    for look in vborrow:
                        if look.stats.station == ii.stats.station:
                            st = Stream(traces=[ii, look])
                            st.merge(fill_value=0)
                            ii = st[0]
                            break
                try:
                    ii.trim(t1, t1 + (process_unit + overlap) * 60,nearest_sample=False, pad=True, fill_value=0)
                    ii.detrend()
                except NotImplementedError:
                    ii.data=np.zeros(points)

                vraw.append(ii)


            for i in h1:
                ii = copy.deepcopy(i)
                if int(float(time_begin['hour'])+runh)%24 == 23:
                    for look in h1borrow:
                        if look.stats.station == ii.stats.station:
                            st = Stream(traces=[ii, look])
                            st.merge(fill_value=0)
                            ii = st[0]
                            break
                try:
                    ii.trim(t1, t1 + (process_unit + overlap) * 60,nearest_sample=False, pad=True, fill_value=0)
                    ii.detrend()
                except NotImplementedError:
                    ii.data = np.zeros(points)

                h1raw.append(ii)

            for i in h2:
                ii = copy.deepcopy(i)
                if int(float(time_begin['hour'])+runh)%24 == 23:
                    for look in h2borrow:
                        if look.stats.station == ii.stats.station:
                            st = Stream(traces=[ii, look])
                            st.merge(fill_value=0)
                            ii = st[0]
                            break
                try:
                    ii.trim(t1, t1 + (process_unit + overlap) * 60,nearest_sample=False, pad=True, fill_value=0)
                    ii.detrend()
                except NotImplementedError:
                    ii.data = np.zeros(points)

                h2raw.append(ii)

            '''
            vwood = []
            h1wood = []
            h2wood = []

            for inv, i, j, k in zip(xmlfiles, h1raw, h2raw, vraw):
                my_tr = Stream(traces=[i, j, k])
                print(k)
                try:
                    my_tr.remove_response(inventory=inv, output="VEL")
                    #my_tr.filter('bandpass', freqmin=0.5, freqmax=10, zerophase=True)
                    my_tr.detrend(type='demean')

                    delta = 1 / sr

                    for i in range(0, len(my_tr)):
                        my_tr[i] = wa(my_tr[i], delta)

                    h1wood.append(my_tr[0].data)
                    h2wood.append(my_tr[1].data)
                    vwood.append(my_tr[2].data)

                except (KeyboardInterrupt, SystemExit):
                    raise

                except:
                    print('fail')
                    h1wood.append(np.ones(points))
                    h2wood.append(np.ones(points))
                    vwood.append(np.ones(points))
            '''

            [h1wood, h2wood, vwood] = parallel_process_traces(xmlfiles, h1raw, h2raw, vraw, sr, points, processes)

            pickle.dump([vfilter, h1filter, h2filter], open(outputdir + day + str(int(float(time_begin['hour'])+runh)%24) + 'wave.p', 'wb'))
            pickle.dump([vwood, h1wood, h2wood], open(outputdir + day + str(int(float(time_begin['hour'])+runh)%24) + 'wood.p', 'wb'))

            vn = []
            hn = []

            count=0
            for i in vfilter:
                a=abs(i)
                if math.isnan(np.median(a)) == True or np.median(a) == 0:
                    a = np.ones(points)
                if len(a) < points:
                    a = np.ones(points)
                if len(a) > points:
                    a=a[0: points]

                for j in range(0,int(round((process_unit+overlap)/1))):
                    if np.median(a[j * 60 * sr:(j + 1) * 60 * sr]) < 10**(-6):
                        a[j * 60 * sr:(j + 1) * 60 * sr] = a[j * 60 * sr:(j + 1) * 60 * sr] ** 0
                    else:
                        #a[j *60*4* sr:(j + 1)*60*4*sr]=a[j*60*4*sr:(j+1)*60*4*sr]/np.median(a[j*60*4*sr:(j+1)*60*4*sr])
                        a[j * 60 * sr:(j + 1) * 60 * sr] = a[j * 60 * sr:(j + 1) * 60 * sr] / np.median(a[j * 60 * sr:(j + 1) * 60 * sr])

                for j in range(0,int(round(points/sr))-60):
                    if np.median(a[j:j + 1 * 60 * sr]) < 10**(-6):
                        a[j*sr:(j+1)*sr] = a[j*sr:(j+1)*sr] ** 0

                a=a**(1/root)
                vn.append(a)

            for i,j in zip(h1filter, h2filter):

                if math.isnan(np.median(i)) == True or np.median(i) == 0:
                    i = np.ones(points)
                if len(i) < points:
                    i = np.ones(points)
                if len(i) > points:
                    i=i[0: points]

                if math.isnan(np.median(j)) == True or np.median(j) == 0:
                    j = np.ones(points)
                if len(j) < points:
                    j = np.ones(points)
                if len(j) > points:
                    j=j[0: points]

                a=(i**2+j**2)**0.5

                for j in range(0,int(round((process_unit+overlap)/1))):
                    if np.median(a[j * 60 * sr:(j + 1) * 60 * sr]) < 10**(-6):
                        a[j * 60 * sr:(j + 1) * 60 * sr] = a[j * 60 * sr:(j + 1) * 60 * sr] ** 0
                    else:
                        #a[j *60*4* sr:(j + 1)*60*4*sr]=a[j*60*4*sr:(j+1)*60*4*sr]/np.median(a[j*60*4*sr:(j+1)*60*4*sr])
                        a[j * 60 * sr:(j + 1) * 60 * sr] = a[j * 60 * sr:(j + 1) * 60 * sr] / np.median(a[j * 60 * sr:(j + 1) * 60 * sr])

                for j in range(0,int(round(points/sr))-60):
                    if np.median(a[j:j + 1 * 60 * sr]) < 10**(-6):
                        a[j*sr:(j+1)*sr] = a[j*sr:(j+1)*sr] ** 0

                a=a**(1/root)
                hn.append(a)

            brp = ssa.med_scan_mul(ptraveltimes, np.array(vn), sr, win, step, processes)
            brmapp = ((brp / (win * sr)) ** 0.5 / nsta) ** root

            brs = ssa.med_scan_mul(straveltimes, np.array(hn), sr, win, step, processes)
            brmaps = ((brs / (win * sr)) ** 0.5 / nsta) ** root

            shape_s = np.shape(brmaps)
            brmapp = brmapp[:, 0:shape_s[1]]
            brmap = np.multiply(brmaps, brmapp)
            brmap = brmap.astype('float32')

            #pickle.dump([brmap], open(outputdir+day+str(int(float(time_begin['hour'])+runh)%24) + 'ssabrmap.p', 'wb'))

            latnum, lonnum = 104, 104
            depnum = 1
            dm = 120
            win = 6
            timestep = 7200

            br3d = np.zeros([timestep + dm, latnum, lonnum], dtype=np.float32)

            for e in range(timestep + dm):
                if e%100 == 0:
                    print(e)

                seedepth = 0

                for i in range(0, latnum):
                    for j in range(0, lonnum):
                        br3d[e, i, j] = brmap[seedepth + (lonnum * i + j) * depnum, e]

            savemat(outputdir + day + str(int(float(time_begin['hour'])+runh)%24) + ".mat", {'cropVol' : br3d})


    print(latgrid)
    print(longrid)

if __name__ == '__main__':
    main()
