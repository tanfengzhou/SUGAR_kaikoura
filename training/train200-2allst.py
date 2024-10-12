from obspy import read
import pickle
import numpy as np

stnm = []
stlas = []
stlos = []
stz = []

with open('./syn_wave/stnm_allst', 'r') as f:
    comein = f.readlines()

for i in comein:
    stnm.append(i[:-1])

with open('./syn_wave/stlas_allst' ,'r') as f:
    comein = f.readlines()

for i in comein:
    stlas.append(float(i))

with open('./syn_wave/stlos_allst', 'r') as f:
    comein = f.readlines()

for i in comein:
    stlos.append(float(i))

with open('./syn_wave/stz_allst', 'r') as f:
    comein = f.readlines()

for i in comein:
    stz.append(float(i))

nsta = len(stnm)
sr = 40
dl = 4 * 60 * sr

readdir = 'syn_wave/EQ1000allst/'
savedir = 'syn_wave/toAI/data1970/02/04clean/'

wf = np.zeros([nsta, 3, (20*60+5)*60*sr])

#ot = np.random.random(4000) * 3600 * 20
#pickle.dump([ot], open(savedir + 'ot.p', 'wb'))

with open(savedir + 'ot.p', 'rb') as f:
    comein = pickle.load(f)

ot = comein[0]

for i in range(4000):
    print(i)
    for j in range(nsta):
        z = read(readdir + 'EQ' + str(i%1000) + '/' + stnm[j] + '.Z.SAC')
        n = read(readdir + 'EQ' + str(i%1000) + '/' + stnm[j] + '.N.SAC')
        e = read(readdir + 'EQ' + str(i%1000) + '/' + stnm[j] + '.E.SAC')
        t_start = int(round(ot[i] * sr))
        t_end = t_start + 4 * 60 * sr
        wf[j][0][t_start:t_end] = wf[j][0][t_start:t_end] + z[0].data[0:4*60*sr]
        wf[j][1][t_start:t_end] = wf[j][1][t_start:t_end] + n[0].data[0:4*60*sr]
        wf[j][2][t_start:t_end] = wf[j][2][t_start:t_end] + e[0].data[0:4*60*sr]

for i in range(nsta):
    z = read(readdir + 'EQ' + str(i) + '/' + stnm[j] + '.Z.SAC')
    n = read(readdir + 'EQ' + str(i) + '/' + stnm[j] + '.N.SAC')
    e = read(readdir + 'EQ' + str(i) + '/' + stnm[j] + '.E.SAC')
    z[0].data = wf[i][0]
    n[0].data = wf[i][1]
    e[0].data = wf[i][2]
    z[0].stats.starttime = '1970-02-04T00:00:00.0Z'
    n[0].stats.starttime = '1970-02-04T00:00:00.0Z'
    e[0].stats.starttime = '1970-02-04T00:00:00.0Z'
    z.write(savedir + stnm[i] + '.Z.SAC')
    n.write(savedir + stnm[i] + '.N.SAC')
    e.write(savedir + stnm[i] + '.E.SAC')

8