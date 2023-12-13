import pickle
import numpy as np
import geopy.distance
from pyproj import Geod
geodetic = Geod(ellps='WGS84')

ymax = 313
ymin = -291
xmax = 400
xmin = -300
zmax = 85
zmin = -1

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

print(latlon2xy(-43.3, 172.0))        # ymax = 179
print(latlon2xy(-41.1, 175.2))          # ymin = -179
print(latlon2xy(-43.2, 175.2))          # xmax = 248
print(latlon2xy(-41.1, 172.0))        # xmin = -106

cutymax = 179+6
cutymin = -179-6
cutxmax = 248+6
cutxmin = -106-6

stnms = []

with open('stnm_fewer', 'r') as f:
    comein = f.readlines()
for i in comein:
    stnms.append(i[:-1])

for stnm in stnms:
    with open('tts_' + stnm + '.p', 'rb') as f:
        ttp = pickle.load(f)

    xyjianju = 0.2
    zjianju = 1

    ttp_cut = ttp[int(round((cutxmin - xmin) / xyjianju)): int(round((cutxmax - xmin) / xyjianju)), int(round((cutymin - ymin) / xyjianju)): int(round((cutymax - ymin) / xyjianju)), :]

    print(np.size(ttp))
    print(np.size(ttp_cut))

    pickle.dump(ttp_cut, open('ttsforlocate_' + stnm + '.p', 'wb'))

7
