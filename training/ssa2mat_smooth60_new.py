from scipy.io import savemat
import pickle
import numpy as np
import math
import geodis
from scipy.ndimage import gaussian_filter

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx


dirread = 'brvideo_soft/noise57_select/'
outputdir = 'training_data_noise57_soft_select/preprocessedDataset/'

studyarea = [35.5, 35.3, 3.58, 4.39]
xygrid = 4         #km
depgrid = 4        #km
margin = 20

latgrid = xygrid/geodis.geodis([[studyarea[0]-0.5,studyarea[1]],[studyarea[0]+0.5,studyarea[1]]])     #degree
longrid = xygrid/geodis.geodis([[studyarea[0],studyarea[1]-0.5],[studyarea[0],studyarea[1]+0.5]])     #degree

studyarea[0] = studyarea[0] - latgrid * margin
studyarea[1] = studyarea[1] - longrid * margin
studyarea[2] = studyarea[2] + latgrid * margin * 2
studyarea[3] = studyarea[3] + longrid * margin * 2

lats = np.arange(studyarea[0],studyarea[0]+studyarea[2],latgrid)
lons = np.arange(studyarea[1],studyarea[1]+studyarea[3],longrid)


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
    ('day','01'),
    ('hour','10'),
    ('min','00'),
    ('second', '00')
])

latnum = len(lats)
lonnum = len(lons)
depnum = 1
dm = 120
win = 6
timestep = 7200

totalhours = int(time_end['hour']) - int(time_begin['hour'])


with open('syn_wave/data1970/01/' + time_begin['day'] + 'clean/input_catalog.csv', 'r') as f:
    comein = f.readlines()

input_catalog = []
for i in comein:
    detail = i.split(',')
    t = float(detail[0])
    lat = float(detail[1])
    lon = float(detail[2])
    dep = float(detail[3])
    mag = float(detail[4])
    input_catalog.append([t, lat, lon, dep, mag])


labels3d = np.zeros([totalhours * timestep + dm, latnum, lonnum], dtype=np.float32)
br3d = np.zeros([totalhours * timestep + dm, latnum, lonnum], dtype=np.float32)

for i in input_catalog:
    step = int(round((i[0] - win / 2)) / 3600 * timestep)
    if step > totalhours * timestep + dm:
        break
    la = find_nearest(lats, i[1])
    lo = find_nearest(lons, i[2])
    labels3d[step,la,lo] = i[4]

labels3d = gaussian_filter(labels3d, sigma=1.5)
norm = 0.0188135
labels3d = labels3d / norm


n = 0
h = 0

with open(dirread + time_begin['day'] + str(h) + 'ssabrmap.p', 'rb') as f:
    comein = pickle.load(f)
brmap2d = comein[0]
brmap2d = brmap2d[:,0:timestep]

for h in range(1,totalhours-1):

    with open(dirread + time_begin['day'] + str(h) + 'ssabrmap.p', 'rb') as f:
        comein = pickle.load(f)
    brmap = comein[0]
    brmap = brmap[:,0:timestep]
    brmap2d = np.concatenate((brmap2d, brmap), axis=1)

h = totalhours-1
with open(dirread + time_begin['day'] + str(h) + 'ssabrmap.p', 'rb') as f:
    comein = pickle.load(f)
brmap = comein[0]
brmap2d = np.concatenate((brmap2d, brmap), axis=1)


for e in range(timestep * totalhours + dm):
    if e % 10000 == 0:
        print(e)
    br3d[e] = brmap2d[:,e].reshape(latnum, lonnum)


for i in input_catalog:
    if i[0] > totalhours * timestep:
        break

    step = int(round((i[0] - win / 2) / 3600 * timestep))
    la = find_nearest(lats, i[1])
    lo = find_nearest(lons, i[2])

    locationchange = np.round(np.random.normal(scale=0.5, size=[3,8])*10)
    locationchange = locationchange.astype(int)

    for j in range(8):
        center = [step + locationchange[0,j], la + locationchange[1,j], lo + locationchange[2,j]]
        if center[0]-30<0 or center[0]+30>timestep*totalhours or center[1]-30<0 or center[1]+30>latnum or center[2]-30<0 or center[2]+30>lonnum:
            continue

        dataone = br3d[center[0]-30:center[0]+30, center[1]-30:center[1]+30, center[2]-30:center[2]+30]
        labelone = labels3d[center[0]-10:center[0]+10, center[1]-10:center[1]+10, center[2]-10:center[2]+10]
        dataone = dataone.astype(np.float32)
        labelone = labelone.astype(np.float32)

        savemat(outputdir + 'labelsTr/' + time_begin['day'] + str(n) + ".mat", {'labelone' : labelone})
        savemat(outputdir + 'imagesTr/' + time_begin['day'] + str(n) + ".mat", {'dataone' : dataone})

        n = n + 1

ran = int(round(n/6))

ran1 = np.round(np.random.random(ran) * timestep * totalhours - 2000)
ran2 = np.round(np.random.random(ran) * (latnum - 60))
ran3 = np.round(np.random.random(ran) * (lonnum - 60))

ran1 = ran1.astype(int)
ran2 = ran2.astype(int)
ran3 = ran3.astype(int)

for i in range(ran):
    center = [ran1[i]+1000, ran2[i]+30, ran3[i]+30]
    dataone = br3d[center[0] - 30:center[0] + 30, center[1] - 30:center[1] + 30, center[2] - 30:center[2] + 30]
    labelone = labels3d[center[0] - 10:center[0] + 10, center[1] - 10:center[1] + 10, center[2] - 10:center[2] + 10]
    dataone = dataone.astype(np.float32)
    labelone = labelone.astype(np.float32)

    savemat(outputdir + 'labelsTr/' + time_begin['day'] + str(n) + ".mat", {'labelone': labelone})
    savemat(outputdir + 'imagesTr/' + time_begin['day'] + str(n) + ".mat", {'dataone': dataone})

    n = n + 1


print(n)


