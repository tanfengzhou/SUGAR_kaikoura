from obspy import geodetics
import numpy as np
import matplotlib.pyplot as plt
import pickle

outputdir = './random_event_allst/'

mc = 3
N = 1000
#aa = np.random.random(N)
#magnitude = - np.log10(aa) + mc
#pickle.dump(magnitude, open(outputdir + 'magnitude.p', 'wb'))

#plt.hist(magnitude, bins=50)
#plt.yscale('log')

with open(outputdir + 'magnitude.p', 'rb') as f:
    comein = pickle.load(f)
magnitude = comein

for i in range(N):
    magnitude[i] = magnitude[i]+1

sourcet = []
for i in magnitude:
    if i>=3 and i<4:
        sourcet.append(np.around((i-3)*0.3+0.1, decimals=3))
    elif i<5:
        sourcet.append(np.around((i-4)*0.6+0.4, decimals=3))
    else:
        sourcet.append(np.around((i-5)*3+1, decimals=3))

sdr = [[254,80,180],[188,45,45],[228,45,90],[148,45,180],[192,45,90],[235,50,180],[213,60,90],[225,70,135],[240,70,180],[220,70,180],[155,45,0],[226,25,135],[219,38,128]] * 77
del sdr[1000]
'''
with open('../gns20200311.csv', 'r') as f:
    comein = f.readlines()

comein.reverse()

evlas = []
evlos = []
evdeps = []

for i in range(17, 917):
    detail = comein[i].split(',')
    evlas.append(detail[5])
    evlos.append(detail[4])
    evdeps.append(detail[7])

for i in range(100):
    evlas.append(np.random.random() * 1.4 + (-42.9))
    evlos.append(np.random.random() * 2.2 + 172.4)
    evdeps.append(np.random.random() * 30)

pickle.dump([evlas, evlos, evdeps], open(outputdir + 'events.p', 'wb'))
'''
with open(outputdir + 'events.p', 'rb') as f:
    comein = pickle.load(f)
[evlas, evlos, evdeps] = comein

evdeps[161] = str(float(evdeps[161]) - 0.01)
evdeps[164] = str(float(evdeps[164]) - 0.01)
evdeps[337] = str(float(evdeps[337]) - 0.01)
evdeps[757] = str(float(evdeps[757]) - 0.01)
evdeps[784] = str(float(evdeps[784]) - 0.01)


for i in range(900, 1000):
    evlas[i] = str(np.around(evlas[i], decimals=6))
    evlos[i] = str(np.around(evlos[i], decimals=6))
    evdeps[i] = str(np.around(evdeps[i], decimals=6))

stlas = []
stlos = []
stz = []
stnm = []

with open('stlas_allst' ,'r') as f:
    comein = f.readlines()

for i in comein:
    stlas.append(float(i))

with open('stlos_allst', 'r') as f:
    comein = f.readlines()

for i in comein:
    stlos.append(float(i))

with open('stz_allst', 'r') as f:
    comein = f.readlines()

for i in comein:
    stz.append(float(i))

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
for x,y,z in zip(evlos, evlas, evdeps):
    plt.scatter(float(x),float(y), cmap="viridis", c=float(z), vmin=0, vmax=30)

#ax.set_xlim([172, 175])
#ax.set_ylim([-43.2, -41.2])
plt.xlabel('lon', fontsize=12)
plt.ylabel('lat', fontsize=12)
cbar = plt.colorbar()
cbar.set_label('Depth (km)', rotation=270, labelpad=30, fontsize=12)
cbar.ax.tick_params(labelsize=12)
plt.tick_params(axis='both', which='both', labelsize=12)
plt.savefig(outputdir + 'events.png')

for e in range(N):

    dist = []
    azi = []

    for i,j in zip(stlas, stlos):
        [d, a, b] = geodetics.base.gps2dist_azimuth(float(evlas[e]), float(evlos[e]), i, j)
        dist.append(d/1000)
        azi.append(a)

    f = open(outputdir + str(e) + '.csv', 'w')
    for i in range(len(stlas)):
        f.write(str(i) + ',' + str(dist[i]) + ',' + str(azi[i]) + ',' + str(stz[i]) + ',' + evdeps[e] + ',' + str(sdr[e][0]) + ','  + str(sdr[e][1]) + ',' + str(sdr[e][2]) + ',' + str(sourcet[e]) + ',' + str(np.around(magnitude[e], decimals=2)) + '\n')

    f.close()

3