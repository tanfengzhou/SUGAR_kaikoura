import numpy as np
import pickle
import taupz
import csv
import geodis
import obspy
from obspy import read
from pathlib import Path
import distance
import math
import copy
import ssa
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.io import savemat

nametag='noise000004'

area='turall'
studyarea = [35.5, 35.3, 3.58, 4.39]
studydepth=[10,10]
xygrid=4         #km
depgrid=4        #km
margin = 20
processes = 120

latgrid=xygrid/geodis.geodis([[studyarea[0]-0.5,studyarea[1]],[studyarea[0]+0.5,studyarea[1]]])     #degree
longrid=xygrid/geodis.geodis([[studyarea[0],studyarea[1]-0.5],[studyarea[0],studyarea[1]+0.5]])     #degree

studyarea[0] = studyarea[0] - latgrid * margin
studyarea[1] = studyarea[1] - longrid * margin
studyarea[2] = studyarea[2] + latgrid * margin * 2
studyarea[3] = studyarea[3] + longrid * margin * 2

time_begin=dict([
    ('year','1970'),
    ('month','01'),
    ('day','01'),
    ('hour','00'),
    ('min','00'),
    ('second', '00')
])

time_end=dict([
    ('year','1970'),
    ('month','01'),
    ('day','08'),
    ('hour','20'),
    ('min','00'),
    ('second', '00')
])

process_unit = 60
overlap = 5             # min

#datadir = 'syn_wave/data1970/02/' + time_begin['day'] + nametag + '/'
outputdir = 'brvideo/' + nametag + '/'

######################################################################

sr=40      # Hz
scanlowf=5
scanhighf=20
root=3
win=6       # seconds
step=0.5      # seconds
points=int(sr*(process_unit+overlap)*60)

######################################################################

with open('../tableP_iasp91_new.p' , 'rb') as f:
    comein=pickle.load(f)
tableP = comein[0]

with open('../tableS_iasp91_new.p' , 'rb') as f:
    comein=pickle.load(f)
tableS = comein[0]

stations=[]
with open('../station_master_all.csv' ,'r', encoding='utf-8-sig') as f:
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

lats = np.arange(studyarea[0],studyarea[0]+studyarea[2],latgrid)
lons = np.arange(studyarea[1],studyarea[1]+studyarea[3],longrid)
deps = np.arange(studydepth[0],studydepth[1]+depgrid,depgrid)

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

    datadir = 'syn_wave/data1970/01/' + day + nametag + '/'

    stlas=[]
    stlos=[]
    stz=[]

    v=[]
    h1=[]
    h2=[]

    for i in station_use:
        name1 = datadir + i[name_index] + '.Z.SAC'
        n = datadir + i[name_index] + '.N.SAC'
        e = datadir + i[name_index] + '.E.SAC'

        if Path(name1).is_file() == True and Path(n).is_file() == True and Path(e).is_file() == True:
            channelz = read(name1)
            v.append(channelz[0])
            channeln = read(n)
            h1.append(channeln[0])
            channele = read(e)
            h2.append(channele[0])

            stlas.append(float(i[lat_index]))
            stlos.append(float(i[lon_index]))
            stz.append(float(i[elev_index]))

    traveldis=[]
    for i in studygrids:
        a=[]
        for j,k in zip(stlas,stlos):
            a.append(distance.dis(j,k,i[0],i[1]))

        traveldis.append(a)

    ptraveltimes=[]
    straveltimes=[]
    for i in range(0,len(traveldis)):
        if i%100==0:
            print(i)
        a=[]
        b=[]
        for j in range(0, len(traveldis[i])):
            timeP = taupz.taupz(tableP, tableS, studygrids[i][2], traveldis[i][j],'P', stz[j], depgrid=0.3, disgrid=0.003)
            a.append(timeP)
            timeS = taupz.taupz(tableP, tableS, studygrids[i][2], traveldis[i][j],'S', stz[j], depgrid=0.3, disgrid=0.003)
            b.append(timeS)

        ptraveltimes.append(a)
        straveltimes.append(b)

    ptraveltimes = np.array(ptraveltimes)
    straveltimes = np.array(straveltimes)


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

    plt.savefig(outputdir +'stations.pdf')
    pickle.dump([stlas, stlos, stz], open(outputdir + 'stations.p', 'wb'))
    
    #######################################################################

    nsta=len(v)
    '''
    if day == time_end['day']:
        duration = int(float(time_end['hour']))
    elif day == time_begin['day']:
        duration = 24 - int(float(time_begin['hour']))
    else:
        duration = 24

    if day == time_end['day'] and day == time_begin['day']:
        duration = int(float(time_end['hour'])) - int(float(time_begin['hour']))
    '''
    if day == '01':
        duration = 10
    elif day == '02':
        duration = 21
    else:
        duration = 20

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
                name1 = datadir + time_begin['month'] + '/' + dayp + '/' + i[name_index] + '.mseed'
                if Path(name1).is_file() == True:
                    channels = read(name1)
                    zcount = 0
                    ncount = 0
                    ecount = 0
                    for j in channels:
                        if j.stats.channel == 'HNZ' and zcount == 0:
                            vborrow.append(j)
                            zcount = 1
                        if (j.stats.channel == 'HN1' or j.stats.channel == 'HNE') and ncount == 0:
                            h1borrow.append(j)
                            ncount = 1
                        if (j.stats.channel == 'HN2' or j.stats.channel == 'HNN') and ecount == 0:
                            h2borrow.append(j)
                            ecount = 1

                    stlas.append(float(i[lat_index]))
                    stlos.append(float(i[lon_index]))
                    stz.append(float(i[elev_index]))

        t1 = UTCDateTime(time_begin['year']+'-'+time_begin['month']+'-'+str(int(float(time_begin['day'])+np.floor(float(time_begin['hour'])+runh)/24))+'T'+str(int(float(time_begin['hour'])+runh)%24)+':00:00')

        vfilter = []
        h1filter = []
        h2filter = []

        for i in v:
            ii = copy.deepcopy(i)
            if int(float(time_begin['hour'])+runh)%24 == 23:
                for look in vborrow:
                    if look.stats.station == ii.stats.station:
                        ii = ii + look
                        break
            try:
                ii.detrend()
                ii.filter(type='bandpass', freqmin=scanlowf, freqmax=scanhighf, corners=3, zerophase=True)
                ii.trim(t1, t1 + (process_unit + overlap) * 60, nearest_sample=False, pad=True, fill_value=0)
                vfilter.append(ii.data)
            except NotImplementedError:
                vfilter.append(np.ones(points))

        for i in h1:
            ii = copy.deepcopy(i)
            if int(float(time_begin['hour']) + runh) % 24 == 23:
                for look in h1borrow:
                    if look.stats.station == ii.stats.station:
                        ii = ii + look
                        break
            try:
                ii.detrend()
                ii.filter(type='bandpass', freqmin=scanlowf, freqmax=scanhighf, corners=3, zerophase=True)
                ii.trim(t1, t1 + (process_unit + overlap) * 60, nearest_sample=False, pad=True, fill_value=0)
                h1filter.append(ii.data)
            except NotImplementedError:
                h1filter.append(np.ones(points))

        for i in h2:
            ii = copy.deepcopy(i)
            if int(float(time_begin['hour']) + runh) % 24 == 23:
                for look in h2borrow:
                    if look.stats.station == ii.stats.station:
                        ii = ii + look
                        break
            try:
                ii.detrend()
                ii.filter(type='bandpass', freqmin=scanlowf, freqmax=scanhighf, corners=3, zerophase=True)
                ii.trim(t1, t1 + (process_unit + overlap) * 60, nearest_sample=False, pad=True, fill_value=0)
                h2filter.append(ii.data)
            except NotImplementedError:
                h2filter.append(np.ones(points))

        '''
        epi=[173.240845, -42.532986]
        epidis=[]
        for i,j in zip(stlas, stlos):
            epidis.append(geodis.geodis([[epi[1],epi[0]],[i,j]]))
        for i in range(0,len(vfilter)):
            plt.plot(vfilter[i]/max(vfilter[i])*3 + epidis[i])
            #plt.plot(h1filter[i]/max(h1filter[i]) + epidis[i])
        '''

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

            for j in range(0,process_unit+overlap):
                a[j *60* sr:(j + 1)*60*sr]=a[j*60*sr:(j+1)*60*sr]/np.median(a[j*60*sr:(j+1)*60*sr])

            #for j in range(0,int(round(process_unit/4))):
            #    a[j *60*4* sr:(j + 1)*60*4*sr]=a[j*60*4*sr:(j+1)*60*4*sr]/np.median(a[j*60*4*sr:(j+1)*60*4*sr])
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

            for j in range(0,process_unit+overlap):
                a[j *60* sr:(j + 1)*60*sr]=a[j*60*sr:(j+1)*60*sr]/np.median(a[j*60*sr:(j+1)*60*sr])

            #for j in range(0,int(round(process_unit/4))):
            #    a[j *60*4* sr:(j + 1)*60*4*sr]=a[j*60*4*sr:(j+1)*60*4*sr]/np.median(a[j*60*4*sr:(j+1)*60*4*sr])
            a=a**(1/root)
            hn.append(a)

        epi = [36.6088,	36.4057]
        epidis = []
        for i, j in zip(stlas, stlos):
            epidis.append(geodis.geodis([[epi[0], epi[1]], [i, j]]))

        for i in range(0,len(vn)):
            plt.plot(vn[i] + epidis[i])
            #plt.plot(hn[i] + i*6+3)

        brp = ssa.med_scan_mul(ptraveltimes, np.array(vn), sr, win, step, processes)
        brmapp = ((brp / (win * sr)) ** 0.5 / nsta) ** root

        brs = ssa.med_scan_mul(straveltimes, np.array(hn), sr, win, step, processes)
        brmaps = ((brs / (win * sr)) ** 0.5 / nsta) ** root

        #brp = ssa.med_scan(ptraveltimes, np.array(vn), sr, win, step)
        #brmapp = ((brp / (win * sr)) ** 0.5 / nsta) ** root

        #brs = ssa.med_scan(straveltimes, np.array(hn), sr, win, step)
        #brmaps = ((brs / (win * sr)) ** 0.5 / nsta) ** root

        shape_s = np.shape(brmaps)
        brmapp = brmapp[:, 0:shape_s[1]]
        brmap = np.multiply(brmaps, brmapp)
        brmap = brmap.astype('float32')

        pickle.dump([brmap], open(outputdir+day+str(int(float(time_begin['hour'])+runh)%24) + 'ssabrmap.p', 'wb'))

        #############################################################################
        # select save format: .p or .mat
        #############################################################################
        '''
        timestep, dm, latnum, lonnum, depnum = 7200, 120, len(lats), len(lons), 1
        br3d = np.zeros([timestep + dm, latnum, lonnum], dtype=np.float32)

        for e in range(timestep + dm):
            if e % 100 == 0:
                print(e)

            seedepth = 0

            for i in range(0, latnum):
                for j in range(0, lonnum):
                    br3d[e, i, j] = brmap[seedepth + (lonnum * i + j) * depnum, e]

        savemat(outputdir + day+str(int(float(time_begin['hour'])+runh)%24) + '.mat', {'cropVol': br3d})
        '''

4
