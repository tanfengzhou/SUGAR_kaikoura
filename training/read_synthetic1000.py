import os

from obspy import read
from obspy.signal import rotate
from obspy.core.utcdatetime import UTCDateTime
import matplotlib.pyplot as plt
import csv

for eq in range(0, 1000):

    print(eq)

    readdir = 'fk_outputs/EQ' + str(eq) + '/'
    outputdir = 'fk_zne/EQ' + str(eq) + '/'

    os.mkdir(outputdir)

    with open(readdir + 'StaInfo.csv', 'r') as f:
        comein = f.readlines()

    azi = []
    dis = []
    idold = []
    id = []

    for i in range(1, len(comein)):
        detail = comein[i].split(',')
        idold.append(int(detail[0]))
        dis.append(float(detail[1]))
        azi.append(detail[2])
        id.append(str(i-1))
        duration = detail[8]

    stlas = []
    stlos = []
    stz = []
    stnm = []

    with open('../station_master_all copy.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')

        for row in reader:
            if row['area'] == 'turall':
                stlas.append(float(row['lat']))
                stlos.append(float(row['lon']))
                stz.append(float(row['elev']))
                stnm.append(row['sta'])

    nsta = len(stnm)

    start = UTCDateTime(0-60)

    for i in range(nsta):
        r = read(readdir + 'WFraw/Sta' + id[i] + '_Azi' + azi[i] + '/WF' + duration + '/sta' + id[i] + '.r')
        t = read(readdir + 'WFraw/Sta' + id[i] + '_Azi' + azi[i] + '/WF' + duration + '/sta' + id[i] + '.t')
        z = read(readdir + 'WFraw/Sta' + id[i] + '_Azi' + azi[i] + '/WF' + duration + '/sta' + id[i] + '.z')

        #plt.plot(r[0].data/max(abs(r[0].data)))
        #plt.plot(t[0].data/max(abs(t[0].data)) + 1)

        n, e = rotate.rotate_rt_ne(r[0].data, t[0].data, 360 - float(azi[i]))
        r[0].data = n
        t[0].data = e

        r.trim(start, start + 240, pad=True, nearest_sample=False, fill_value=0)
        t.trim(start, start + 240, pad=True, nearest_sample=False, fill_value=0)
        z.trim(start, start + 240, pad=True, nearest_sample=False, fill_value=0)

        #plt.plot(r[0].data / max(abs(r[0].data)) + 2)
        #plt.plot(t[0].data / max(abs(t[0].data)) + 3)

        z.write(outputdir + stnm[idold[i]] + '.Z.SAC')
        r.write(outputdir + stnm[idold[i]] + '.N.SAC')
        t.write(outputdir + stnm[idold[i]] + '.E.SAC')

    4
