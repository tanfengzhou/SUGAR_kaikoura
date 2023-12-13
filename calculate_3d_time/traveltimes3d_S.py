import pickle
import numpy as np
import geopy.distance
from pyproj import Geod
geodetic = Geod(ellps='WGS84')

ymax = 313                    # km
ymin = -291
xmax = 400
xmin = -300
zmax = 85
zmin = -1
xyjianju = 0.2                # km
zjianju = 1                   # km
originlat = -41.7638
originlon = 172.9037

def xy2latlon(x, y):
    y = -y
    angle = np.rad2deg(np.arctan2(y, x))
    angle = 130 - angle
    dis = np.sqrt(x ** 2 + y ** 2)
    c = geopy.distance.distance().destination((originlat, originlon), bearing=angle, distance=dis)
    return([c.latitude, c.longitude])

def latlon2xy(lat, lon):
    az12, az21, dist = geodetic.inv(originlon, originlat, lon, lat)
    angle = 130 - az12
    y = - np.sin(np.deg2rad(angle)) * dist / 1000
    x = np.cos(np.deg2rad(angle)) * dist / 1000
    return([x, y])


stnms = []
stlas = []
stlos = []
stz = []

with open('stnm_fewer', 'r') as f:
    comein = f.readlines()
for i in comein:
    stnms.append(i[:-1])

with open('stlas_fewer', 'r') as f:
    comein = f.readlines()
for i in comein:
    stlas.append(float(i))

with open('stlos_fewer', 'r') as f:
    comein = f.readlines()
for i in comein:
    stlos.append(float(i))

with open('stz_fewer', 'r') as f:
    comein = f.readlines()
for i in comein:
    stz.append(float(i))

with open('studygrids.p', 'rb') as f:
    comein = pickle.load(f)
studygrids = comein[0]

print(len(studygrids))

with open('vs3d.p', 'rb') as f:
    comein = pickle.load(f)
vp3d = comein

ptraveltimes=[]

tt_station = {}

for i in stnms:
    with open('tts_' + i + '.p', 'rb') as f:
        tttable = pickle.load(f)
    tt_station[i] = tttable

n = 0
for i in studygrids:
    a = []
    n = n + 1
    if n % 100 == 0:
        print(n)
    for j,k,l,m in zip(stlas, stlos, stz, stnms):
        [hypox, hypoy] = latlon2xy(i[0], i[1])
        ttp = tt_station[m]
        traveltime = ttp[int(round((hypox - xmin) / xyjianju)), int(round((hypoy - ymin) / xyjianju)), int(round((i[2] - zmin) / zjianju))]
        [stationx, stationy] = latlon2xy(j, k)
        v = vp3d[int(round((stationx - xmin) / xyjianju)), int(round((stationy - ymin) / xyjianju)), 0]
        traveltime = traveltime - l / 1000 / v
        a.append(traveltime)

    ptraveltimes.append(a)

ptraveltimes = np.array(ptraveltimes)
ptraveltimes_save = ptraveltimes.astype('float32')

pickle.dump(ptraveltimes_save, open('straveltimes3d_0.2.p', 'wb'))

4
