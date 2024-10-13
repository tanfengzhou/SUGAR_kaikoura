from obspy import read
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


totalhour = 20
noise_std = np.arange(0.00002, 0.0002, (0.0002-0.00002)/totalhour)
#noise_std = np.ones(totalhour) * 0.0001

z = []
n = []
e = []

readdir = 'syn_wave/toAI/data1970/02/03clean/'
savedir = 'syn_wave/toAI/data1970/02/03/'


for i in range(nsta):
    st = read(readdir + stnm[i] + '.Z.SAC')
    for h in range(totalhour-1):
        nn = noise_std[h]
        st[0].data[60*60*sr*h:60*60*sr*(h+1)] = st[0].data[60*60*sr*h:60*60*sr*(h+1)] + np.random.normal(0, nn, 60*60*sr)

    h = totalhour - 1
    nn = noise_std[h]
    st[0].data[60*60 * sr * h:len(st[0].data)] = st[0].data[60*60 * sr * h:len(st[0].data)] + np.random.normal(0, nn ,60*65 * sr)
    st.write(savedir + stnm[i] + '.Z.SAC')

#for i in range(len(z)):
#    plt.plot(z[i].data/max(abs(z[i].data))+i)

for i in range(nsta):
    st = read(readdir + stnm[i] + '.N.SAC')
    for h in range(totalhour-1):
        nn = noise_std[h]
        st[0].data[60*60*sr*h:60*60*sr*(h+1)] = st[0].data[60*60*sr*h:60*60*sr*(h+1)] + np.random.normal(0, nn, 60*60*sr)

    h = totalhour-1
    nn = noise_std[h]
    st[0].data[60*60 * sr * h:len(st[0].data)] = st[0].data[60*60 * sr * h:len(st[0].data)] + np.random.normal(0, nn ,60*65 * sr)
    st.write(savedir + stnm[i] + '.N.SAC')


for i in range(nsta):
    st = read(readdir + stnm[i] + '.E.SAC')
    for h in range(totalhour-1):
        nn = noise_std[h]
        st[0].data[60*60*sr*h:60*60*sr*(h+1)] = st[0].data[60*60*sr*h:60*60*sr*(h+1)] + np.random.normal(0, nn, 60*60*sr)

    h = totalhour-1
    nn = noise_std[h]
    st[0].data[60*60 * sr * h:len(st[0].data)] = st[0].data[60*60 * sr * h:len(st[0].data)] + np.random.normal(0, nn ,60*65 * sr)
    st.write(savedir + stnm[i] + '.E.SAC')



8